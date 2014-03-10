#!/usr/bin/env perl
# Download an SRA file, dump it to fastq, and shuffle the reads (if PE)
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Temp qw/tempdir/;
use File::Basename;

sub logmsg{print STDERR "@_\n"; }
exit main();

sub main{
  my $settings={
    checkexecs=>1,
  };
  GetOptions($settings,qw(help tempdir=s checkexecs!)) or die $!;
  die usage() if(!@ARGV);
  $$settings{tempdir}||=tempdir( CLEANUP => 1 );
  logmsg "Temporary directory is $$settings{tempdir}";

  my $query=join(" ",@ARGV);

  checkForEdirect() if($$settings{checkexecs});
  my $SRA=findSraId($query,$settings);
  for my $sra(@$SRA){
    my $fastq=downloadSra($sra,$settings);
    $fastq=shuffleReads($fastq,$settings);
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

sub downloadSra{
  my($acc,$settings)=@_;
  logmsg "Downloading accession $acc from SRA";
  system("fastq-dump -I --split-files -O $$settings{tempdir} -v -v -v '$acc' 1>&2");
  die if $?;
  my @fastq=("$$settings{tempdir}/${acc}_1.fastq","$$settings{tempdir}/${acc}_2.fastq");#glob("$$settings{tempdir}/$acc*.fastq");
  logmsg "Created files ".join(" ",@fastq); # TODO: are these sorted properly?
  return \@fastq;
}

sub shuffleReads{
  my($fastq,$settings)=@_;
  my $acc=fileparse($$fastq[0],"_1.fastq");
  my $out="$$settings{tempdir}/$acc.shuffled.fastq";
  logmsg "Shuffling the reads into $out";
  system("run_assembly_shuffleReads.pl ".join(" ",@$fastq)." > $out");
  die if $?;
  return $out;
}

sub printReads{
  my($fastq,$settings)=@_;
  logmsg "Outputting the fastq to stdout";
  system("cat '$fastq'"); die if $?;
  return 1;
}

sub usage{
  "Download an SRA file, dump it to fastq, and shuffle the reads.  Expects PE Fastq
  Warning: only one Illumina run is expected. Multiple runs have not been tested.
  Usage: $0 query text > out.fastq
  -t tempdir Default: one will be made for you
  -nocheck to not check for executables
  "
}
