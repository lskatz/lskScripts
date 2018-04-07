#!/usr/bin/env perl

# Convert a fastq to a fasta/qual combo using BioPerl, with some Linux commands

use Bio::Perl;
use Data::Dumper;
use strict;
use warnings;
use threads;
use Thread::Queue;
use Getopt::Long;

my $settings={};

$|=1;
my %numSequences; # static for a subroutine

exit(main());

sub main{
  my $usage="Usage: $0 -i inputFastqFile [-n numCpus -q outputQualfile -f outputFastaFile]";
  die($usage) if(@ARGV<1);

  GetOptions($settings,('numCpus=s','input=s','qualOut=s','fastaOut=s','help'));
  die $usage if($$settings{help});

  my $file=$$settings{input}||die("input parameter missing");
  my $outfasta=$$settings{fastaOut}||"$file.fasta";
  my $outqual=$$settings{qualOut}||"$file.qual";
  my $numCpus=$$settings{numCpus}||1;

  my @subfile=splitFastq($file,$numCpus);
  for my $f(@subfile){
    threads->create(\&convert,$f,"$f.fasta","$f.qual");
  }
  $_->join for (threads->list);

  # join the sub files together
  joinFastqFiles(\@subfile,$file);

  return 1;
}

sub convert{
  my($file,$outfasta,$outqual)=@_;

  my $numSequences=numSequences($file);
  my $reportEvery=int($numSequences/100) || 1;
  print "$numSequences sequences to convert in $file\n";

  my $in=Bio::SeqIO->new(-file=>$file,-format=>"fastq-illumina");
  my $seqOut=Bio::SeqIO->new(-file=>">$outfasta",-format=>"fasta");
  my $qualOut=Bio::SeqIO->new(-file=>">$outqual",-format=>"qual");
  my $seqCount=0;
  my $percentDone=0;
  while(my $seq=$in->next_seq){
    $seqOut->write_seq($seq);
    $qualOut->write_seq($seq);
    $seqCount++;
    if($seqCount%$reportEvery == 0){
      $percentDone++;
      print "$percentDone%..";
    }
  }
  print "Done with subfile $file.\n";
  return 1;
}

sub joinFastqFiles{
  my($subfile,$outfileBasename)=@_;
  my($command,$subfasta,$subqual);

  # fasta
  $subfasta.="$_.fasta " for(@$subfile);
  $command="cat $subfasta > $outfileBasename.fasta";
  system($command);
  
  # qual
  $subqual.="$_.qual " for (@$subfile);
  $command="cat $subqual > $outfileBasename.qual";
  system($command);

  return 1;
}

sub splitFastq{
  my($file,$numCpus)=@_;
  my $prefix="FQ"; # for fastq
  my $numSequences=numSequences($file);
  my $numSequencesPerFile=int($numSequences/$numCpus);
  my $numSequencesPerFileRemainder=$numSequences % $numCpus;
  my $numLinesPerFile=$numSequencesPerFile*4; # four lines per read
  system("rm -r tmp;mkdir tmp;");
  system("split -l $numLinesPerFile $file 'tmp/FQ'");

  return glob("tmp/FQ*");
}


sub numSequences{
  my $file=shift;
  return $numSequences{$file} if($numSequences{$file});
  my $num=`grep -c '^\@' $file`;
  chomp($num);
  $numSequences{$file}=$num;
  return $num;
}
