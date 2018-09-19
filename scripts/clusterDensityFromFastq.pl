#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use Getopt::Long qw/GetOptions/;
use POSIX qw/ceil/;
use threads;
use Thread::Queue;

local $0 = basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i tile-size|size-of-tile=s)) or die $!;
  $$settings{'tile-size'} ||= 1;
  $$settings{numcpus}||=1;

  my $fastqPerThread = ceil(scalar(@ARGV) / $$settings{numcpus});

  my $printQ = Thread::Queue->new();
  $printQ->enqueue("File\tspots-per-mm2");

  my @thr;
  for(my $i=0;$i<$$settings{numcpus};$i++){
    my @fastq = splice(@ARGV,0,$fastqPerThread);
    $thr[$i] = threads->new(sub{
      my($fastqArr, $printQ)=@_;
      for my $fastq(@$fastqArr){
        logmsg $fastq;
        my $density = clusterDensity($fastq,$settings);
        $printQ->enqueue([$fastq, $density]);
      }
      return scalar(@fastq);
    }, \@fastq, $printQ);
  }

  # start the printer
  my $printerThread = threads->new(\&printer, $printQ);

  # join the threads
  for(@thr){
    $_->join;
  }

  # Terminate multithreaded printing
  $printQ->enqueue(undef);
  $printerThread->join();

  return 0;
}

sub printer{
  my($Q)=@_;
  while(defined(my $toPrint = $Q->dequeue)){
    if(ref($toPrint) eq 'ARRAY'){
      print join("\t",@$toPrint)."\n";
    } else {
      print "$toPrint\n";
    }
  }
}

sub clusterDensity{
  my($fastq,$settings)=@_;

  # Try to get some compile-time speedup
  my $colonRegex = qr/:/;
  my $whitespaceRegex = qr/\s+/;

  my %tileCount;
  open(my $fh, "-|", "zcat $fastq") or die "ERROR: could not zcat $fastq: $!";
  while(my $header=<$fh>){
    # Burn three lines. We're only looking at the header.
    <$fh>;
    <$fh>;
    <$fh>;

    chomp($header);
    my($firstPart,undef) = split($whitespaceRegex, $header);
    my($instrument, $run, $flowcell, $lane, $tile, $x, $y) = split($colonRegex, $firstPart);
    #my(undef, undef, undef, undef, $tile, $x, $y) = split($colonRegex, $firstPart);
    $tileCount{$tile}++;
  }
  close $fh;

  my $total=0;
  while(my($tile,$count)=each(%tileCount)){
    $total+=$count;
  }
  my $averagePerTile = sprintf("%0.2f", $total/scalar(keys(%tileCount))/$$settings{'tile-size'});

  return $averagePerTile;
}

sub usage{
  "$0: Calculates cluster density of a fastq file with casava style headers
  Usage: $0 [options] *.fastq.gz
  
  --tile-size  1  Size of the tile in square mm.
  --numcpus    1
  "
}
