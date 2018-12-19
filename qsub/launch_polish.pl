#!/usr/bin/env perl
#$ -pe smp 1-16
#$ -S /bin/bash
#$ -cwd -V
#$ -o polish.log -j y
#$ -N polishAssembly

use 5.12.0;
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw/GetOptions/;
use File::Basename qw/basename/;
use Env::Modulecmd;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use File::Copy qw/cp/;

Env::Modulecmd::purge();
Env::Modulecmd::load(qw(nanopolish/0.8.3 minimap2/2.7 Python/2.7.13 
  samtools/1.4.1 minimap2/2.7 bowtie2/2.3.3.1 pilon/1.22));

my $numcpus=$ENV{NSLOTS} || 24;

exit(main());

sub main{

  my $settings={};
  GetOptions($settings,qw(help nanopore=s illumina1=s illumina2=s outfile|output|outfasta=s)) or die $!;
  $$settings{nanopore}//="";
  $$settings{illumina1}//="";
  $$settings{illumina2}//="";
  $$settings{outfile}//="polished.fasta";
  $$settings{tempdir}//=tempdir("asm_polish_XXXXXX", TMPDIR=>1, CLEANUP=>1);

  my $ref=shift @ARGV;
  die usage() if(!$ref || $$settings{help});

  if($$settings{nanopore}){
    $ref=nanopolish($ref,$$settings{nanopore},$settings);
  }
  if($$settings{illumina1}){
    $ref=pilon     ($ref,$$settings{illumina1},$$settings{illumina2},$settings);
  }

  # Easy to print results with cat!
  system("cat $ref");
  die if $?;
  
  return 0;
}

sub nanopolish{
  ...;
}

sub pilon{
  my($ref, $R1, $R2, $settings)=@_;

  my $symlinkRef="$$settings{tempdir}/unpolished.fasta";
  cp($ref, $symlinkRef);
  $ref=$symlinkRef;

  # People say to run Pilon four times in a row, single threaded
  my $maxRuns=4;
  for(my $i=1;$i<=4;$i++){
    my $bam="$$settings{tempdir}/pilon$i.sorted.bam";

    system("bowtie2-build $ref $ref"); die if $?;
    system("bowtie2 -x $ref -1 $R1 -2 $R2 -p $numcpus | samtools sort -T $$settings{tempdir}/samtoolssort -o $bam");
    die if $?;
    system("samtools index $bam"); die if $?;

    system("pilon --genome '$ref' --frags $bam --output $$settings{tempdir}/pilon$i --changes --threads $numcpus --fix snps,indels");
    die if $?;

    # update for the next iteration
    $ref="$$settings{tempdir}/pilon$i.fasta";
  }

  return $ref;
}

sub usage{
  "Polish an assembly with nanopore or illumina reads.\n".
  basename($0)." [--nanopore=nanopore.fastq] --illumina1=R1.fastq --illumina2=R2.fastq --outfile=out.fasta contigs.fasta
  "
}

