#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use Getopt::Long qw/GetOptions/;
use Bio::Kmer;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i kmerlength=s)) or die $!;
  $$settings{kmerlength}||="29,17";
  $$settings{numcpus}||=1;
  die usage() if($$settings{help} || @ARGV < 2);

  my @kmerlength = sort {$a <=> $b} split(/,/,$$settings{kmerlength});
  for(my $i=0;$i<@ARGV;$i+=2){
    my $R1=$ARGV[$i];
    my $R2=$ARGV[$i+1];

    my $entries = readFastq($R1,$R2,$settings);
    my $topKmers = countKmers($R1, \@kmerlength, $settings);
  }

  return 0;
}

sub readFastq{
  my($R1,$R2,$settings)=@_;

  my %entry;

  open(my $fh1,"gzip -cd $R1|") or die "ERROR: could not open $R1: $!";
  open(my $fh2,"gzip -cd $R2|") or die "ERROR: could not open $R2: $!";
  while(my $id1 = <$fh1>){
    my $seq1  = <$fh1>;
    my $plus1 = <$fh1>;
    my $qual1 = <$fh1>;
    my $id2   = <$fh2>;
    my $seq2  = <$fh2>;
    my $plus2 = <$fh2>;
    my $qual2 = <$fh2>;
    my $combinedSeq = $seq1."N".$seq2;
    $combinedSeq=~s/\s+//g;

    chomp $id1;
    $entry{$id1} = {
      entry => "$id1\n$seq1$plus1$qual1$id2$seq2$plus2$qual2",
      seq   => $combinedSeq,
    }
  }
  return \%entry;
}

sub countKmers{
  my($fastq, $kmerLength, $settings) = @_;
  
  for my $k(@$kmerLength){
    my $kmerObj = Bio::Kmer->new($fastq, {kmerlength=>$k, numcpus=>$$settings{numcpus}, gt=>2});
    my $kmers = $kmerObj->kmers;
    die Dumper $kmers;
  }
}

sub usage{
  local $0 = basename $0;
  "Sorts fastq entries by common kmers of varying lengths
  Usage: $0 R1.fastq.gz R2.fastq.gz
  --kmerlength  29,17
  --numcpus     1
  "
}

