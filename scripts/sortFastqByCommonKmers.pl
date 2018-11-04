#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use Getopt::Long qw/GetOptions/;
#use List::MoreUtils qw/uniq/;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help kmerlength=s)) or die $!;
  $$settings{kmerlength}||="29,17";
  die usage() if($$settings{help} || @ARGV < 2);

  my @kmerlength=split(/,/,$$settings{kmerlength});
  for(my $i=0;$i<@ARGV;$i+=2){
    my $R1=$ARGV[$i];
    my $R2=$ARGV[$i+1];

    my $kmerInfo      = countKmers ($R1,$R2,\@kmerlength,$settings);
    my $sortedEntries = sortEntries($kmerInfo,$settings);
  }

  return 0;
}

sub sortEntries{
  my($kmerInfo, $settings)=@_;
  my $entries     = $$kmerInfo{read};
  my $kmerCounter = $$kmerInfo{kmerCounter};
  my $kmer_to_read= $$kmerInfo{kmer_to_read};

  # gzip seems to do better if kmer lengths are ascending
  my @kmerLength = sort {$a <=> $b} keys(%$kmerCounter);

  # Find the most common kmer of each length
  my %topKmer;
  for my $kmerLength(@kmerLength){
    my $maxcount=0;
    my $topKmer ="";
    while(my($kmer, $count) = each(%{$$kmerCounter{$kmerLength}})){
      if($count > $maxcount){
        $topKmer  = $kmer;
        $maxcount = $count;
      }
    }
    $topKmer{$kmerLength} = $topKmer;
  }

  # For each kmer length, print out resulting reads
  for my $kmerLength(@kmerLength){
    my @readId = @{ $$kmer_to_read{$topKmer{$kmerLength}} };

    # print the reads
    for my $id(@readId){
      if($$entries{$id}){
        print $$entries{$id};
        # remove from the entries hash
        $$entries{$id} = 0;
      }
    }
  }

  # Print the rest
  for my $entry(sort {$a cmp $b} values(%$entries)){
    print $entry;
  }
}

sub countKmers{
  my($R1,$R2,$kmerlength,$settings)=@_;
  my @kmerlength=@$kmerlength; # bring into scope to make it faster(?)

  # read{id} => 8-line-str, 
  my %read;
  # kmerCounter => {31 =>{AAA=>4,...}, 21=>{...}
  my %kmerCounter;
  # kmer_to_read => {AAA => id1, ...}
  my %kmer_to_read;

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
    my $entry = "$id1$seq1$plus1$qual1$id2$seq2$plus2$qual2";
    chomp($id1,$seq1,$plus1,$qual1,$id2,$seq2,$plus2,$qual2);

    # Save the read
    $read{$id1} = $entry;

    # kmerizing
    for my $k(@kmerlength){
      my $seqkLength = length($seq1) - $k + 1;
      
      # kmers for each window
      for(my $i=0;$i<$seqkLength;$i++){
        my $kmer = substr($seq1, $i, $k);

        # count the kmer
        $kmerCounter{$k}{$kmer}++;

        # Mark that this kmer is in this read
        #$kmer_to_read{$kmer}{$id1}=1;
        push(@{$kmer_to_read{$kmer}}, $id1);
        #$kmer_to_read{$kmer} = [uniq(@{$kmer_to_read{$kmer}})];
      }
    }
  }
  close $fh2;
  close $fh1;
  
  return {read=>\%read,kmerCounter=>\%kmerCounter,kmer_to_read=>\%kmer_to_read};
}


sub usage{
  local $0 = basename $0;
  "Sorts fastq entries by common kmers of varying lengths
  Usage: $0 R1.fastq.gz R2.fastq.gz
  --kmerlength  29,17
  "
}

