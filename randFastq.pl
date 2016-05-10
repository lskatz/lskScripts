#!/usr/bin/env perl
# Randomizes the order of fastq reads
#

use strict;
use warnings;
use Getopt::Long qw/GetOptions/;
use Data::Dumper qw/Dumper/;
use File::Basename qw/fileparse/;
use List::Util qw/shuffle/;

local $0=fileparse $0;
sub logmsg { print "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help pe|paired-end freq|frequency=f)) or die $!;
  $$settings{freq}||=1;
  die usage() if($$settings{help});

  my @fastq=@ARGV;
  die usage() if(!@fastq);

  my $reads=readFastqs(\@fastq,$settings);

  printRandomReads($reads,$settings);

  return 0;
}

sub readFastqs{
  my($fastq,$settings)=@_;

  my $linesPerEntry=4;
  if($$settings{pe}){
    $linesPerEntry=8;
  }

  # Get this out of the hash in case it helps with speed
  my $freq=$$settings{freq};

  my @reads;
  for my $f(@$fastq){
    my($name,$dir,$ext)=fileparse($f,qw(.gz));
    my $fastqFh;
    if($ext eq '.gz'){
      open($fastqFh,"zcat $f |") or die "ERROR: could not zcat $f for reading: $!";
    } else {
      open($fastqFh,$f) or die "ERROR: could not open $f: $!";
    }
    while(my $entry=<$fastqFh>){
      for(2..$linesPerEntry){
        $entry.=<$fastqFh>;
      }

      # Randomly skip reads if a random number is greater
      # than the user-defined threshold.
      next if(rand() > $freq);

      push(@reads,$entry);
    }
    close $fastqFh;
  }

  return \@reads;
}

sub printRandomReads{
  my($reads,$settings)=@_;
  
  for my $entry(shuffle(@$reads)){
    print $entry;
  }
}

sub usage{
  "$0: randomize the order of reads in a fastq file
  Usage: $0 file.fastq[.gz] [file2.fastq...] > rand.fastq

  --paired-end      If the file is interleaved
  --frequency   1   Frequency of reads to keep (values: 0-1)
  "
}

