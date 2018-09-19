#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use Getopt::Long qw/GetOptions/;

local $0 = basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help tile-size|size-of-tile=s)) or die $!;
  $$settings{'tile-size'} ||= 1;

  print "File\tspots-per-mm2\n";
  for my $fastq(@ARGV){
    logmsg $fastq;
    my $density = clusterDensity($fastq,$settings);
    print "$fastq\t$density\n";
  }

  return 0;
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
  "
}
