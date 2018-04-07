#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Perl;
use Bio::Restriction::EnzymeCollection;
use Bio::Restriction::Analysis;
use Data::Dumper;
use File::Basename qw/basename/;
use Getopt::Long;

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n"}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help outtype=s enz|enzyme=s@)) or die $!;
  $$settings{enz}||=[qw(AscI)];
  $$settings{outtype}||="bed";

  my @seq=@ARGV;
  die usage() if($$settings{help}||!@seq);

  my $all_collection = Bio::Restriction::EnzymeCollection->new(); 
  my @enz;
  for my $enzName(@{ $$settings{enz} }){
    my $re=$all_collection->get_enzyme($enzName); 
    # Try capitalizing the enzyme or other tricks
    # if it isn't found right away.
    if(!$re){
      my $EnzName=ucfirst($enzName);
      $re=$all_collection->get_enzyme($EnzName);
      logmsg "Tried transforming $enzName to $EnzName";
    }
    if(!$re){
      my $enzN4me=$enzName;
      $enzN4me=~s/1(I*)$/I$1/g;  # replace tail ones with Is
      $re=$all_collection->get_enzyme($enzN4me);
      logmsg "Tried transforming $enzName to $enzN4me";
    }
    if(!$re){
      die "ERROR: I do not understand enzyme $enzName";
    }
    push(@enz,$re);
  }

  my $counter=0;
  for my $gbk(@seq){
    my $in=Bio::SeqIO->new(-file=>$gbk); 
    while(my $seq=$in->next_seq){ 
      my $seqLength=$seq->length;
      my $ra=Bio::Restriction::Analysis->new(-seq=>$seq);

      if($$settings{outtype} eq 'sizes'){
        my @fragments=map{length($$_{seq})} $ra->fragment_maps(@enz);
        for(@fragments){
          print join("\t",$seq->id,$_)."\n";
        }
      } elsif($$settings{outtype} eq 'bed'){
        for my $re(@enz){
          my @pos=$ra->positions($re->name);
          for my $pos(@pos){
            # I'm not 100% sure why this position math works, but
            # it matches up with what Apollo genome browser does.
            my $start=$pos-$re->cut+1;
            my $end=$start+$re->recognition_length-1;
            print join("\t",$seq->id,$start,$end,$re->name.++$counter)."\n";
          }
        }
      } else {
        die "ERROR: I don't understand outtype $$settings{outtype}";
      }
    }
  }

  return 0;
}

sub usage{
  local $0=basename $0;
  "Usage: $0 *.fasta > restrictionAnalysis.bed
  --enzyme    AscI  The enzyme to digest with. Can suppy
                    multiple --enzyme arguments.
  --outtype   bed   Outputs a bed file of cut size coordinates.
                    If 'sizes' is supplied instead, then 
                    fragment sizes will be output.
  "
}

