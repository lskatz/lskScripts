#!/usr/bin/env perl
# Find the minimum abundance of kmers
# Original script was in Python at 
#   https://gist.github.com/alexjironkin/4ed43412878723491240814a0d5a6ed6/223dea45d70c9136703a4afaab0178cdbfbd2042
# Original author of python script: @alexjironkin, Public Health Englad
# I wanted to increase capatibility with Perl and have a more
# standalone script instead of relying on the khmer package.
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename fileparse/;
use List::Util qw/max/;

use constant MAXINT =>  ~0;
use constant MININT => -~0;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help kmerlength|kmer=i delta=i gt|greaterthan=i hist|histogram=s)) or die $!;
  $$settings{kmerlength}||=21;
  $$settings{delta}     ||=100;
  $$settings{gt}        ||=0;

  my($fastq)=@ARGV;
  die usage() if(!$fastq || $$settings{help});
  die "ERROR: I could not find fastq at $fastq" if(!-e $fastq);

  print join("\t",qw(kmer count))."\n";

  my $histogram=[];
  if($$settings{hist} && -e $$settings{hist}){
    $histogram=readHistogram($$settings{hist},$settings);
  } else {
    my $kmercount=countKmers($fastq,$$settings{kmerlength},$settings);
    $histogram=kmerToHist($kmercount,$settings);
  }
  my $valley   =findTheValley($histogram,$$settings{delta},$settings);
  print join("\t",@$valley)."\n";

  return 0;
}

sub readHistogram{
  my($infile,$settings)=@_;
  my @hist=(0);
  open(HIST,$infile) or die "ERROR: could not read $infile: $!";
  logmsg "Reading histogram from $infile";
  while(<HIST>){
    chomp;
    my($count,$countOfCounts)=split /\t/;
    $hist[$count]=$countOfCounts;
  }
  close HIST;

  # ensure defined values
  for(my $i=0;$i<@hist;$i++){
    $hist[$i] //= 0;
  }

  return \@hist;
}

sub countKmers{
  my($fastq,$kmerlength,$settings)=@_;
  my %kmer=();

  # pure perl to make this standalone
  #open(FASTQ,"<", $fastq) or die "ERROR: could not read $fastq: $!";
  open(FASTQ,"zcat $fastq | ") or die "ERROR: could not read $fastq: $!";
  my $i=0;
  while(my $id=<FASTQ>){
    my $seq=<FASTQ>;
    chomp($seq);
    my $numKmersInRead=length($seq)-$kmerlength+1;

    # Count kmers in a sliding window.
    # We must keep this loop optimized for speed.
    for(my $j=0;$j<$numKmersInRead;$j++){
      $kmer{substr($seq,$j,$kmerlength)}++;
    }

    # Burn the quality score lines
    <FASTQ>;
    <FASTQ>;

  }
  close FASTQ;

  return \%kmer;
}

sub kmerToHist{
  my($kmercountHash,$settings)=@_;
  my %hist=();
  my @hist=(0);

  for my $kmercount(values(%$kmercountHash)){
    $hist{$kmercount}++;
  }

  # Turn this hash into an array
  for(1..max(keys(%hist))){
    $hist[$_] = $hist{$_} || 0;
  }

  if($$settings{hist}){
    logmsg "Writing histogram to $$settings{hist}";
    open(HIST, ">", $$settings{hist}) or die "ERROR: could not write histogram to $$settings{hist}: $!";
    for(my $i=0;$i<@hist;$i++){
      print HIST "$i\t$hist[$i]\n";
    }
    close HIST;
  }

  return \@hist;
}

sub findTheValley{
  my($hist, $delta, $settings)=@_;

  my($min,$max)=(MAXINT,MININT);
  my($minPos,$maxPos)=(0,0);
  my @maxTab=();
  my @minTab=();

  my $lookForMax=1;

  my $numZeros=0; # If we see too many counts of zero, then exit.
  
  for(my $kmerCount=$$settings{gt}+1;$kmerCount<@$hist;$kmerCount++){
    my $countOfCounts=$$hist[$kmerCount];
    if($countOfCounts == 0){ 
      $numZeros++;
    }
    if($countOfCounts > $max){
      $max=$countOfCounts;
      $maxPos=$kmerCount;
    }
    if($countOfCounts < $min){
      $min=$countOfCounts;
      $minPos=$kmerCount;
    }

    if($lookForMax){
      if($countOfCounts < $max - $delta){
        push(@maxTab,[$maxPos,$max]);
        $min=$countOfCounts;
        $minPos=$kmerCount;
        $lookForMax=0;
      }
    }
    else{
      if($countOfCounts > $min + $delta){
        push(@minTab,[$minPos,$min]);
        $max=$countOfCounts;
        $maxPos=$kmerCount;
      }
    }

    last if($numZeros > 3);
  }

  return $minTab[0];
}


sub usage{
  "$0: Find the valley between two peaks on a set of kmers
       such that you can discard the kmers that are 
       likely representative of contamination.
  Usage: $0 file.fastq.gz
  --gt     0   Look for the first peak at this kmer count
               and then the next valley.
  --kmer   21  kmer length
  --delta  100 How different the counts have to be to
               detect a valley or peak
  --hist   ''  A file to write the histogram, or a file
               to read the histogram if it already exists.
               Useful if you want to rerun this script.
  "
}

