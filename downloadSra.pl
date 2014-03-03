#!/usr/bin/env perl
# Download an SRA file, dump it to fastq, and shuffle the reads (if PE)
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Temp qw/tempdir/;

sub logmsg{print STDERR "@_\n"; }
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help pairedEnd tempdir=s)) or die $!;
  die usage() if(!@ARGV);
  $$settings{tempdir}||=tempdir( CLEANUP => 1 );
  logmsg "Temporary directory is $$settings{tempdir}";

  my $query=join(" ",@ARGV);

  checkForEdirect();
  my $SRA=findSraId($query,$settings);
  for my $sra(@$SRA){
    my $fastq=downloadSra($sra,$settings);
    $fastq=shuffleReads($fastq,$settings) if($$settings{pairedEnd});
    printReads($fastq,$settings);
  }
  return 0;
}

sub checkForEdirect{
  for my $exec(qw(esearch efetch xtract fastq-dump run_assembly_shuffleReads.pl cat)){
    system("which $exec >& /dev/null");
    die "ERROR: could not find $exec in your PATH" if $?;
  }
  return 1;
}

# (for i in `seq 68 77`; do CFSAN="CFSAN100$i"; esearch -db sra -query $CFSAN|efetch -format docsum|xtract -element Runs; done;) | perl -lane '/acc="(.+?)"/; $acc=$1; print "$acc";'
sub findSraId{
  my($query,$settings)=@_;
  my $xml=`esearch -db sra -query '$query' | efetch -format docsum | xtract -element Runs`;
  die if $?;

  my @acc;
  while($xml=~/acc="(.+?)"/g){
    push(@acc,$1);
  }
  die "Could not find any accession pertaining to the query! Query was\n  $query" if(!@acc);
  logmsg "Found the following accessions: ".join(", ",@acc);

  return \@acc;    
}

#system("fastq-dump -I --split-files -v -v $acc"); die "Problem getting $CFSAN - $acc" if $?
sub downloadSra{
  my($acc,$settings)=@_;
  system("fastq-dump -I --split-files -O $$settings{tempdir} -v -v -v '$acc'");
  die if $?;
  my @fastq=glob("$$settings{tempdir}/$acc*.fastq");
  logmsg "Created files ".join(" ",@fastq); # TODO: are these sorted properly?
  return \@fastq;
}

sub shuffleReads{
  my($fastq,$settings)=@_;
  ...
}

sub printReads{
  my($fastq,$settings)=@_;
  ...
}

sub usage{
  "Download an SRA file, dump it to fastq, and shuffle the reads (if PE)
  Warning: only one Illumina run is expected. Multiple runs have not been tested.
  Usage: $0 query text > out.fastq
  -p to indicate that you expect paired-end
  "
}
