#!/usr/bin/env perl

use warnings;
use strict;
use Bio::SeqIO;
use Bio::FeatureIO;
use Bio::Restriction::Enzyme;
use Bio::Restriction::EnzymeCollection;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use FindBin;

my $enzCollection=Bio::Restriction::EnzymeCollection->new();

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help force-type=s fragments tempdir=s enzyme=s@ out=s)) or die $!;
  $$settings{tempdir}||=tempdir("addCutSites.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{enzyme}//=[];
  my $richseq=$ARGV[0];

  die usage() if(!$richseq || $$settings{help});

  logmsg "Temporary directory is $$settings{tempdir}";
  my $featHash=getCutSiteFeatures($richseq,$settings);

  my $richSeqHash=addFeaturesToRichseq($richseq,$featHash,$settings); 
  addFragments($richSeqHash,$settings) if($$settings{fragments});

  my $outFh;
  if($$settings{out}){
    $outFh=Bio::SeqIO->new(-file=>$$settings{out});
  } else {
    $outFh=Bio::SeqIO->new(-fh=>\*STDOUT,-format=>"genbank");
  }
  for my $seq(values(%$richSeqHash)){
    $outFh->write_seq($seq);
  }

  return 0;
}

sub getCutSiteFeatures{
  my($richseq,$settings)=@_;

  my $enzArg="";
  for(@{$$settings{enzyme}}){
    $enzArg.="--enzyme $_ ";
  }
  $enzArg=~s/^\s+|\s+$//g; # whitespace trim

  logmsg "Running pfge on the genome";
  system("perl $FindBin::RealBin/pfgeOnGenome.pl $enzArg $richseq > $$settings{tempdir}/cuts.bed");
  die "ERROR with pfgeOnGenome.pl" if $?;

  logmsg "Reading cut sites into memory";
  my $featIn=Bio::FeatureIO->new(-file=>"$$settings{tempdir}/cuts.bed");
  my %feat;
  while(my $feat=$featIn->next_feature){
    my $name=($feat->get_tag_values("Name"))[0];
    $feat->display_name($name); 
    # Allowed tags:
    # http://song.cvs.sourceforge.net/viewvc/song/ontology/so.obo
    # What kind of feature is it?
    # See if the Bio module can tell us type I, II, or III
    my $type;
    my $rec_length=eval{
      my $enzName=$name;
      $enzName=~s/\d+$//;
      my $enzObj=$enzCollection->get_enzyme($enzName);
      $type=$enzObj->type;
      return $enzObj->recognition_length;
    } || "";
    if($@){
      logmsg "WARNING: I either can't get the type, recognition site, or the actual enzyme from the name $name: $@";
    }
    
    # Give it a good primary tag name
    my $primary_tag;
    if($$settings{'force-type'}){
      $primary_tag=$$settings{'force-type'};
    } elsif($rec_length==6){
      $primary_tag="six_cutter_restriction_site";
    } elsif($rec_length==8){
      $primary_tag="eight_cutter_restriction_site";
    } elsif($rec_length==4){
      $primary_tag="four_cutter_restriction_site";
    } elsif($type){
      $primary_tag="type_${type}_enzyme_restriction_site";
    } else{
      $primary_tag="restriction_enzyme_recognition_site"; # most generic
    }
    $feat->primary_tag($primary_tag);  

    my $seqid=$feat->seq_id; # save this value before removing
    $feat->remove_tag($_) for(qw(source seq_id)); # ridiculous, but these are set for printing
    push(@{$feat{$seqid}},$feat);
  }

  return \%feat;
}

sub addFragments{
  my($richSeqHash,$settings)=@_;

  logmsg "Adding fragments in-memory";

  my $enzymeRegex=join('|',@{$$settings{enzyme}});
  #my $enzymeRegex='\Q'.join('\E|\Q',@{$$settings{enzyme}}).'\E';

  my $fragmentCounter=0;
  while(my($seqid,$seq)=each(%$richSeqHash)){
    my $startPos=1;
    my @feat=$seq->get_SeqFeatures;
    for my $f(@feat){
      next if($f->name!~/$enzymeRegex/);

      my $endPos=$f->start;
      my $fragLength=($endPos-$startPos+1);
      my $name="Fragment_$fragLength";
      my $fragment=Bio::SeqFeature::Generic->new(
        -start          => $startPos,
        -end            => $endPos,
      );
      $fragment->add_tag_value('note',"Length: $fragLength");
      $fragment->add_tag_value('Name',$name);
      $fragment->add_tag_value('locus_tag',$name);
      if($$settings{'force-type'}){
        $fragment->primary_tag($$settings{'force-type'});
      } else {
        $fragment->primary_tag("restriction_fragment");
      }
      $seq->add_SeqFeature($fragment);

      $startPos=$endPos+1;
    }
    # Add the last fragment from the last cut site to the end of the seq
    my $fragLength=($seq->length-$startPos+1);
    my $name="Fragment_$fragLength";
    my $endFragment=Bio::SeqFeature::Generic->new(
      -start          => $startPos,
      -end            => $seq->length,
    );
    $endFragment->add_tag_value('note',"Length: $fragLength");
    $endFragment->add_tag_value('Name',$name);
    $endFragment->add_tag_value('locus_tag',$name);
    if($$settings{'force-type'}){
      $endFragment->primary_tag($$settings{'force-type'});
    } else {
      $endFragment->primary_tag("restriction_fragment");
    }
    $seq->add_SeqFeature($endFragment);
  }
}

sub addFeaturesToRichseq{
  my($richseq,$featHash,$settings)=@_;
  
  logmsg "Adding new cut site features into sequences in memory";
  my %seqOut;
  my $in=Bio::SeqIO->new(-file=>$richseq);
  while(my $seq=$in->next_seq){
    $$featHash{$seq->id}//=[];
    for(@{$$featHash{$seq->id}}){
      $seq->add_SeqFeature($_); 
    }
    $seqOut{$seq->id}=$seq;
  }
  return \%seqOut;
}

sub usage{
  local $0=basename $0;
  "$0: Adds restriction cut sites to a richseq file
  Usage: $0 --enzyme AscI in.gbk > out.gbk
  --enzyme     AscI    A restriction enzyme to annotate.
                       Multiple --enzyme allowed.
  --fragments          Add fragment tracks
  --out        ''      If an output file is given, then
                       write to that file. Otherwise,
                       will write to stdout.
  --force-type ''      Force cut site features to be
                       whatever you choose, eg, CDS
  "
}
