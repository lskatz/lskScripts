#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/basename/;
use Data::Dumper;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  
  my($infile)=@ARGV;
  die usage() if(!$infile || $$settings{help});
  
  my @genomeList=`cut -f 5 $infile | sort| grep . | uniq`;
  die if $?;
  chomp(@genomeList);
  SNPs_allToVcf($infile,\@genomeList,$settings);
  return 0;
} 

sub SNPs_allToVcf{
  my($infile,$genomeList,$settings)=@_;

  my @nt=qw(A C G T);
  my $numGenomes=scalar(@$genomeList);

  local $0=basename $0;
  print "##fileformat=VCFv4.1\n##source=kSNP3, $0\n";
  print join("\t",'##CHROM',qw(POS ID REF ALT QUAL FILTER INFO FORMAT),@$genomeList)."\n";

  # Get a list of zeros for the GTs later on
  my @zeroVariantTags=(0) x $numGenomes; #split(//, "0" x $numGenomes);

  # index the genomeList
  my %genomeIndex;
  @genomeIndex{@$genomeList}=keys(@$genomeList);

  my @GT=@zeroVariantTags;
  my %altIndex=();
  my @ALT=('.');
  my($id,$kmer,$variant,$x,$genome)=('.') x 5;
  open(IN,$infile) or die "ERROR: could not read $infile for reading: $!";
  while(<IN>){
    s/^\s+|\s+$//g; # whitespace trim
    
    # Whenever there is a blank line, the ID is about to increment.
    # Print the VCF line.
    if(/^$/){
      #print Dumper \@ALT;
      #next if(join("",@GT) eq join("",@zeroVariantTags)); # don't print non-variants
      next if(scalar(@ALT) < 2);
      next if( (grep{$_ > 0} @GT) < 2);

      # Print the VCF line.
      shift(@ALT); # remove the dot.
      print join("\t",'.','.',$kmer,'.',join(",",@ALT), '.', 'PASS', "NS=$numGenomes", 'GT', @GT)."\n";

      # reset
      @GT=@zeroVariantTags;
      %altIndex=();
      @ALT=('.');
      next;
    }
    
    ($id,$kmer,$variant,$x,$genome)=split(/\t/,$_);
    $id=uc($id);
    $variant=uc($variant);

    if(!defined($altIndex{$variant})){
      $altIndex{$variant}=scalar(@ALT);
      push(@ALT,$variant);
    }
    $GT[$genomeIndex{$genome}]=$altIndex{$variant};
  }
  close IN;
}

sub usage{
  local $0=basename $0;
  "$0: transform a kSNP3 output into a vcf file
  Usage: $0 kSNP3.out/SNPs_all > kSNP3.vcf
  SNPs_all file is formatted with the tab-delimited fields
  ID  kmer  variant  x  genomeName
  "
}
