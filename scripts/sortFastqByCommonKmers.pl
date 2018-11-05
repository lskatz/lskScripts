#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use Getopt::Long qw/GetOptions/;
use List::MoreUtils qw/uniq/;
use Bio::Kmer;

local $0 = basename $0;

sub logmsg{print STDERR "$0: @_\n"}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i kmerlength=s mincount=i)) or die $!;
  $$settings{kmerlength}||="51,29,17";
  $$settings{numcpus}||=1;
  $$settings{mincount}||=10;
  die usage() if($$settings{help} || @ARGV < 2);

  my @kmerlength = sort {$a <=> $b} split(/,/,$$settings{kmerlength});
  for(my $i=0;$i<@ARGV;$i+=2){
    my $R1=$ARGV[$i];
    my $R2=$ARGV[$i+1];

    logmsg "Reading fastq files";
    my $entries  = readFastq($R1,$R2,$settings);
    logmsg "Finding common kmers in R1";
    my $topKmers = topKmers($R1, \@kmerlength, $settings);
    logmsg "Finding common kmers in R2";
    push(@$topKmers, 
                @{ topKmers($R2, \@kmerlength, $settings) }
    );
    $topKmers = [uniq(@$topKmers)];
    logmsg "Using ".scalar(@$topKmers)." kmers to sort";
    logmsg "Printing reads by sorted kmer cluster";
    sortAndPrintReads($entries, $topKmers, $settings);
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
    my $combinedSeq = $seq1."~".$seq2;
    $combinedSeq=~s/\s+//g;

    chomp $id1;
    $entry{$id1} = {
      entry => "$id1\n$seq1$plus1$qual1$id2$seq2$plus2$qual2",
      seq   => $combinedSeq,
    };
  }
  return \%entry;
}

sub topKmers{
  my($fastq, $kmerLength, $settings) = @_;
  my $minCount = $$settings{mincount} || 10;
  
  my @topKmer;
  for my $k(@$kmerLength){
    logmsg "Counting $k-mers";
    my $kmerObj = Bio::Kmer->new($fastq, {kmerlength=>$k, numcpus=>$$settings{numcpus}, gt=>2});
    my $kmers = $kmerObj->kmers;
    while(my($kmer,$count) = each(%$kmers)){
      next if($count < $minCount);
      push(@topKmer,$kmer);
      logmsg "Saving kmer length $k with count $count: $kmer";
    }
  }

  return \@topKmer;
}

sub sortAndPrintReads{
  my($entries, $topKmers, $settings)=@_;

  my @sortedKmers = sort{
    length($a) cmp length($b)
    || $a cmp $b
  } @$topKmers;

  for my $kmer (@sortedKmers){
    logmsg "Printing kmers associated with $kmer";
    my @sortedId = sort{$$entries{$a} cmp $$entries{$b}} keys(%$entries);
    for my $id(@sortedId){
      # Use index() as a faster alternative to regex
      if(index($$entries{$id}{seq},$kmer) != -1){
        # Print the entry and then remove its listing
        print $$entries{$id}{entry};
        delete($$entries{$id});
      }
    }
  }

  # Print the rest
  my @sortedId = sort{$$entries{$a} cmp $$entries{$b}} keys(%$entries);
  for my $id(@sortedId){
    print $$entries{$id}{entry};
  }

  return 1;
}

sub usage{
  "Sorts fastq entries by common kmers of varying lengths
  Usage: $0 R1.fastq.gz R2.fastq.gz
  --kmerlength  29,17
  --numcpus     1
  "
}

