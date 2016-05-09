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
  GetOptions($settings,qw(help pe|paired-end)) or die $!;
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
  Usage: $0 file.fastq [file2.fastq...] > rand.fastq

  --paired-end      If the file is interleaved
  "
}

