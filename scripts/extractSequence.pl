#!/usr/bin/env perl

# Author: Lee Katz <lkatz@cdc.gov>
# Debugging help from Taylor Griswold <TGriswold@cdc.gov>

use strict;
use warnings;
use Bio::Perl;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use List::Util qw/min max/;
use Bio::SeqFeature::Generic;

my $settings={
  contig=>[],
};

exit(main());

sub main{
  die usage() if(@ARGV<1);
  GetOptions($settings,('infile=s', 'start=s', 'nameOfOrganism=s', 'end=s', 'contig=s@', 'outfile=s', 'help'));
  die usage() if($$settings{help});

  my $file=$$settings{infile} || die "Missing infile\n".usage();
  my $end||=$$settings{end};
  my $outfile=$$settings{outfile} or die "Missing outfile\n".usage();
  my $start=$$settings{start} || 1;
  $start=1 if($start<1);

  my $seq=extractSeqs($file,$start,$settings);
  
  my $seqout=Bio::SeqIO->new(-file=>">$outfile");
  $seqout->write_seq($_) for(@$seq);
  
  print "Output file is in $outfile\n";

  return 1;
}

sub extractSeqs{
  my($file,$start,$settings)=@_;
  my @return;

  # get the correct contig Seq
  my $targetSeq=[];
  my $revcom=[];
  if($$settings{contig}){
    my %seqHash=seqHash($file);
    for my $contigId(@{ $$settings{contig} }){
      if($contigId=~/^\-(.+)$/){
        $contigId=$1;
        push(@$revcom,1);
      } else {push(@$revcom,0);}

      die "Could not find contig $contigId in $file\n" if(!$seqHash{$contigId});
      push(@$targetSeq,$seqHash{$contigId});
    }
  } else {
    my $seqin=Bio::SeqIO->new(-file=>$file);
    $targetSeq=[$seqin->next_seq];
    $$settings{contig}=[map($_->id,@$targetSeq)];
  }

  # actually perform extraction
  for(my $i=0;$i<@$targetSeq;$i++){
    my $seq=$$targetSeq[$i];
    my $rc=$$revcom[$i];
    my $tmp=extractSeq($seq,$start,$rc,$settings);
    push(@return,$tmp);
  }
  return \@return;
}

sub extractSeq{

  my($targetSeq,$start,$revcom,$settings)=@_;
  my $end=$$settings{end}||$targetSeq->length;

  my @feat=$targetSeq->get_SeqFeatures();

  $start=1 if($start<1);
  $end=$targetSeq->length if($end>$targetSeq->length);

  print "Extracting and highlighting features in contig ".$targetSeq->id." from $start to $end\n";
  my $truncSeq=$targetSeq->trunc($start,$end);  
  $truncSeq=$truncSeq->revcom if($revcom);
  
  # removes all previously determined features from truncSeq because of potential duplicate features
  $truncSeq->remove_SeqFeatures;
  
  my $seqLength=$end-$start+1;
  my $organism=($targetSeq->species)?
    join(" ",$targetSeq->species->binomial):
    join(" ","StrainXYZ");
  $organism=$$settings{nameOfOrganism} if($$settings{nameOfOrganism});
  # add in the source feature
  $truncSeq->add_SeqFeature(
    Bio::SeqFeature::Generic->new(
      -start=>1,
      -end=>$seqLength,
      -primary=>'source',
      -tag=>{
        organism=>$organism,
      },
    )
  );

  # sort the features appropriately
  @feat=sort {
    return $b->location->start <=> $a->location->start if $revcom;
    return $a->location->start <=> $b->location->start;
  } @feat;

  for my $feat(@feat){
	
    my $subseqFeatStart=$feat->location->start-$start+1;
    my $subseqFeatEnd=$feat->location->end-$start+1;
    my $strand=$feat->strand;

    if($revcom){
      $subseqFeatStart=($end-$feat->location->start+1);
      $subseqFeatEnd=($end-$feat->location->end+1);
      ($subseqFeatStart,$subseqFeatEnd)=($subseqFeatEnd,$subseqFeatStart);
      $strand=$strand*-1;
    } else {
     
    }
    next if( ($subseqFeatStart<1 && $subseqFeatEnd<1) || ($subseqFeatStart>$seqLength && $subseqFeatEnd>$seqLength) );

    my @gene=($feat->has_tag("gene"))?$feat->get_tag_values('gene'):();
    print join("_",$subseqFeatStart,$subseqFeatEnd,$strand,@gene)."\n";

    next if(ref($feat->location) ne "Bio::Location::Simple");
	
    $subseqFeatStart=1 if($subseqFeatStart<1);
    $subseqFeatEnd=$seqLength if($subseqFeatEnd>$seqLength);


    $feat->location->start($subseqFeatStart);
    $feat->location->end($subseqFeatEnd);
    $feat->strand($strand);
    $truncSeq->add_SeqFeature($feat);
  }
  
  #$truncSeq->seq(lc($truncSeq->seq)); # lc the sequence so that it can be highlighted
  #highlightFeatures($truncSeq,$settings);

  return $truncSeq;
}

sub highlightFeatures{
  my($seq,$settings)=@_;
  my $highlight=$$settings{highlight_type}||"uc";
  # TODO return false if not a rich seq
  for my $feat($seq->get_SeqFeatures){
    next if($feat->primary_tag eq 'source');
    my($start,$end)=($feat->location->start,$feat->location->end);
    if($highlight eq 'uc'){
      my $firstPartOfTheSequence=($start>1)?$seq->subseq(1,$start-1):"";
      my  $lastPartOfTheSequence=($end<$seq->length)?$seq->subseq($end+1,$seq->length):"";
      my $sequence=$firstPartOfTheSequence.uc($seq->subseq($start,$end)).$lastPartOfTheSequence;
      $seq->seq($sequence);
    }
  }
}

# makes a hash of sequences, given a fasta file
sub seqHash{
  my $filename=shift;
  my($seqio,%seq);
  if(-e $filename){
    $seqio=Bio::SeqIO->new(-file=>$filename);
  }
  # sequence object
  elsif(ref($filename)=~/seq/i) {
    die ref($filename)." is not supported yet";
  }
  else{
    die "Cannot determine the type of sequence for $filename";
  }
  while(my $seq=$seqio->next_seq){
    my $id=$seq->id;
    $seq{$id}=$seq;

    # make a second key without genome info and zero padding
    $id=~s/^.+?_0*//;
    $seq{$id}=$seq;
  }
  return %seq;
}


sub usage{
  "Usage: perl $0 -i inputGenomicFile -s start -e end
  -i input file: the file extension dictates the format
  -o outfile
  -s start coordinate. 1-based
  -e end coordinate. 1-based
  -c contig 
    put a negative sign in front of a contig to indicate revcom
    you may have multiple -c args
  -n name of organism (useful for genome browsers)
  ";
}
