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
  GetOptions($settings,qw(help tempdir=s checkexecs! only-paired-end only-print-SRRs!)) or die $!;
  die usage() if(!@ARGV || $$settings{help});
  $$settings{'only-print-SRRs'}||=0;
  $$settings{tempdir}||=tempdir( CLEANUP => 1 );
  logmsg "Temporary directory is $$settings{tempdir}";

  my $query=join(" ",@ARGV);

  checkForEdirect() if($$settings{checkexecs});
  my $SRA=findSraId($query,$settings);
  if($$settings{'only-print-SRRs'}){
    print join("\t",$query,@$SRA)."\n";
    return 0;
  }
  ACCESSION: 
  for my $sra(@$SRA){
    my $fastq=downloadSra($sra,$settings);
    
    # If the user only wants paired end but this is single end,
    # don't let anything be printed. Move on.
    if($$settings{"only-paired-end"}){
      for(@$fastq){
        if(!-f $_){
          logmsg "WARNING: user only wanted paired end; however, only single end was detected for accession $sra\n  I will not print reads from accession $sra.";
          next ACCESSION;
        }
      }
    }

    $fastq=shuffleReads($fastq,$settings);
    printReads($fastq,$settings);

    # remove any traces
    #unlink $_ for(glob("$$settings{tempdir}/*.fastq"));
  }
  return 0;
}

sub checkForEdirect{
  for my $exec(qw(esearch efetch xtract fastq-dump prefetch run_assembly_shuffleReads.pl cat)){
    system("which $exec 2>/dev/null");
    die "ERROR: could not find $exec in your PATH" if $?;
  }
  return 1;
}

sub findSraId{
  my($query,$settings)=@_;
  my $numTries=0;
  my $command="esearch -db sra -query '$query' | efetch -format docsum | xtract -element Runs";
  my $xml=`$command`;
  while($? && $numTries++ < 20){
    logmsg "Command failed! $!\n Command was\n   $command";
    $xml=`$command`;
    sleep int(rand(4)) + 1; # sleep for a random time in case other scripts are colliding
  }
  die "ERROR: tried 20 times and failed! $!" if $?;

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
  my $numTries=0;
  logmsg "Downloading accession $acc from SRA";
  system("prefetch '$acc' 1>&2");
  # Do not check for errors on prefetch because fastq-dump will work with or without it

  logmsg "Converting to fastq format";
  # "fastq-dump --defline-seq '$seqIdTemplate' --defline-qual '+' --split-files -O $tmpdir --gzip $F{srarun_acc} "
  $$settings{seqIdTemplate}||='@$ac_$sn[_$rn]/$ri';
  my $command="fastq-dump --defline-seq '$$settings{seqIdTemplate}' --defline-qual '+' --split-files -O $$settings{tempdir}  '$acc' 1>&2";
  system($command);
  while($? && $numTries++ < 20){
    logmsg "Command failed! $!\n Command was\n   $command";
    system($command);
    sleep int(rand(4)) + 1; # sleep for a random time in case other scripts are colliding
  }
  die "ERROR: tried fastq-dump 20 times and failed! $!" if $?;

  my @fastq=("$$settings{tempdir}/${acc}_1.fastq","$$settings{tempdir}/${acc}_2.fastq");#glob("$$settings{tempdir}/$acc*.fastq");
  logmsg "Created files ".join(" ",@fastq); # TODO: are these sorted properly?
  return \@fastq;
}

sub shuffleReads{
  my($fastq,$settings)=@_;
  my $acc=fileparse($$fastq[0],"_1.fastq");
  system("ls -lhS $$settings{tempdir} 1>&2");

  # see if the output is paired end or not.
  # If not, then just return the single end.
  if(!-f "$$settings{tempdir}/$acc"."_2.fastq"){
    logmsg "WARNING: ".$acc."_2.fastq was not found. Returning single-end read set";
    return $$fastq[0];
  }

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
  --only-paired-end to not accept single end runs.
  --only-print-SRRs Do not download, but do print the accessions
  "
}
