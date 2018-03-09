#!/usr/bin/env perl

# splits a genbank file from gendb back out into individual loci (ie contigs or chromosomes)
# Author: Lee Katz <lkatz@cdc.gov>
# With a major assist from Alex Luo <luo.chengwei@gatech.edu>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use Bio::Perl;

my $thisFile=fileparse($0);
exit(main());

sub main{
  my $settings={
    linkerSeq=>"NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN",
  };
  die usage($settings) if(@ARGV<1);
  GetOptions($settings,qw(input=s output=s linkerSeq=s debug help));
  die usage($settings) if($$settings{help});
  for my $param(qw(input output linkerSeq)){
    die "Error: need param $param\n".usage($settings) if(!$$settings{$param});
  }

  splitGbk($$settings{input},$$settings{output},$settings);

  print "Gbk is in $$settings{output}\n";
  return 0;
}

sub splitGbk{
  my($in,$out,$settings)=@_;
  my $seq=Bio::SeqIO->new(-file=>$in)->next_seq;
  my $outseq=Bio::SeqIO->new(-file=>">$out");
  my @feat=$seq->get_SeqFeatures;
  my $desc=$seq->desc;

  # filter for only these features
  # TODO add in all features to accept. The danger is in adding the source feature, which should be split between contigs. I can worry about that another day.
  @feat=grep({$_->primary_tag=~/gene|cds|tRNA|rRNA/i} @feat);

  # split the gbk by the linker sequence
  my $linkerLength=length($$settings{linkerSeq});
  my $linker="\Q$$settings{linkerSeq}\E";
  my @seq=split(/$linker/i,$seq->seq);
  
  print( $seq->length-(scalar(@seq-1)*$linkerLength)." is the expected length\n") if($$settings{debug});

  # write the sequences
  my $totalLength=0;
  my $contigBoundary=1;
  my $contigStart=1;
  my $i=0; # contig number
  
  foreach my $s (@seq){
    $i++;

    # don't ignore that there is a linker sequence
    # add linker: 1-L L-2-L ... L-(n-1)-L  L-n
    my $contigLength=length($s);
    if($i==1){                                            # first contig
      $s.=$$settings{linkerSeq};
    }elsif($i==($#seq)){                                  # last contig
      $s=$$settings{linkerSeq}.$s.$$settings{linkerSeq};
    }else{                                                # middle contigs
      $s=$$settings{linkerSeq}.$s.$$settings{linkerSeq};
    }
    
    my $seqObj=Bio::Seq->new(-id=>$i,-seq=>$s,-desc=>$desc);
    my $seqLength=$contigLength;
    $contigBoundary=length($s);
    print "$contigStart-$contigBoundary\n" if($$settings{debug});

    # alter gene locations according to the contig boundaries
    while(@feat>0){
      my $f=shift(@feat);
    # look at the current feature(s). Are they in between contigs?
      # TODO split the gene between contigs
      if($f->location->start < 1){
        next;
      }

      # if location out of range, unshift it back so that it is on the next contig
      my $j=0;
      $j=1 if ($i==1);
      
      my $genePosOffset=$totalLength+($i+$j-2)*$linkerLength;
      if(
          ($f->location->end-$genePosOffset) > $contigBoundary 
        #&& $f->location->start-$genePosOffset>=$contigBoundary
        ){
        unshift(@feat,$f);
        last;
      }
      
      $f->location->start($f->location->start-$genePosOffset);
      $f->location->end($f->location->end-$genePosOffset);
      
      $seqObj->add_SeqFeature($f);
    }
    
    $totalLength+=$seqLength;

    # add seq features
    $outseq->write_seq($seqObj);
  }
  print "Final derived length: $totalLength\n" if($$settings{debug});
  
  return 0;
}

sub usage{
  my $settings=shift;
  "Splits a genbank file from gendb back out into individual loci
  Usage: $thisFile -i in -o out
  -i input genbank file
  -o output genbank file
    required: .gbk extension
  -d
    debug mode
  -l linker
    default: $$settings{linkerSeq}
  -h
    this help menu
  ";
}
