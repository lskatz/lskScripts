#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use Getopt::Long;
use threads;
use Thread::Queue;

local $0=basename $0;
sub logmsg {print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i)) or die $!;
  $$settings{numcpus}||=1;

  die usage() if(!@ARGV || $$settings{help});
  
  # Find the counts for each fastq file
  my %counts;

  my $Q=Thread::Queue->new(@ARGV);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $Q->enqueue(undef);
    $thr[$_]=threads->new(\&getQualityCounts,$Q,$settings);
  }

  logmsg "All fastq files enqueued!";
  for my $t(@thr){
    logmsg "Waiting to join thread ".$t->tid;
    my $c = $t->join();
    while(my($nt, $positionalQualityCounts)=each(%$c)){
      while(my($pos, $qualityCounts)=each(%$positionalQualityCounts)){
        while(my($quality, $count)=each(%$qualityCounts)){
          $counts{$nt}{$pos}{$quality}+=$count;
        }
      }
    }
  }
  #print Dumper \%counts;die;

  # Print an ART format
  my @nt=sort(keys(%counts));
  for my $nt (@nt){
    my $positionalQualityCounts=$counts{$nt};
    my @posArr=sort {$a<=>$b} keys(%$positionalQualityCounts);

    for my $pos (@posArr){
      my @qualArr=sort {$a cmp $b} keys(%{$$positionalQualityCounts{$pos}});

      # Print the first row for this position, the quality value
      print "$nt\t$pos";
      for my $quality (@qualArr){
        #print "\t$quality";
        print "\t";
        print ord($quality)-33;
      }
      print "\n";

      # Print the second row for this position, the counts at the quality value
      print "$nt\t$pos";
      for my $quality (@qualArr){
        print "\t".$$positionalQualityCounts{$pos}{$quality};
      }
      print "\n";
    }

  }
  
  return 0;
}

sub getQualityCounts{
  my($Q,$settings)=@_;

  my %counts;
  while(defined(my $fastq=$Q->dequeue)){

    logmsg "Processing $fastq";

    my $lineCounter=0;
    open(my $fh, "zcat $fastq | ") or die "ERROR: could not read $fastq: $!";
    while(<$fh>){
      my $seq=<$fh>;
      chomp($seq);
      my @seq  = split(//, $seq);
      <$fh>; # burn the plus line
      my @qual = split(//, scalar(<$fh>));
      # I _could_ chomp @qual but the last item won't be
      # reached if I do the for loop right.

      my $readLength=@seq;
      for(my $i=0;$i<$readLength;$i++){
        $counts{$seq[$i]}{$i}{$qual[$i]}++;
      }
    }
    close $fh;
  }

  return \%counts;
}

sub usage{
  "Usage: $0 *_R1.fastq.gz > profile_R1.txt
          $0 *_R2.fastq.gz > profile_R2.txt
  --numcpus   1   Number of cpus
  "
}
