#!/usr/bin/env perl
# Uses a tree to reads directory and a Snp-Pipeline run to determine
# true positives, false negatives, false positives
#

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help snppipeline=s ttr=s));
  die usage() if($$settings{help});
  die "ERROR: need Snp-Pipeline project name\n".usage() if(!$$settings{snppipeline});
  die "ERROR: need TreeToReads project name\n".usage() if(!$$settings{ttr});

  # where are the Snp-Pipeline SNPs?
  my $snppipelineSnps=snppipelineSnps($$settings{snppipeline},$settings);

  # where are the real SNPs?
  my $realSnps=ttrSnps($$settings{ttr},$settings);
  
  # compare
  my($TP,$TN,$FP,$FN)=compareSnps($snppipelineSnps,$realSnps,$settings);

  # Report
  print join("\t",qw(S-P TTR TP TN FP FN))."\n";
  print join("\t",$$settings{snppipeline},$$settings{ttr},$TP, $TN, $FP, $FN)."\n";
}

sub snppipelineSnps{
  my($proj,$settings)=@_;
  my %pos;

  my $list="$proj/snplist.txt";
  open(SNPMATRIX,"<",$list) or die "ERROR: could not open $list for reading: $!";
  while(<SNPMATRIX>){
    next if(/^#/);
    chomp;
    my($chr,$pos,$count,@genomes)=split /\t/;
    $pos{$pos}=1;
  }
  close SNPMATRIX;

  return \%pos;
}

sub ttrSnps{
  my($proj,$settings)=@_;
  my %pos;

  my $matrix="$proj/var_site_matrix";
  open(TTRMATRIX,"<",$matrix) or die "ERROR: could not open $matrix for reading: $!";
  while(<TTRMATRIX>){
    chomp;
    my($chr,$ref,$pos)=split /\s+/;
    $pos{$pos}=$ref;
  }
  close TTRMATRIX;

  return \%pos;
}

sub compareSnps{
  my ($snppipelineSnps,$realSnps,$settings)=@_;

  # Initialize counts to zero
  my($TP,$TN,$FP,$FN)=split(//,"0" x 4);
  
  # How many of the real SNPs were found?
  while(my($truePos,$trueRef)=each(%$realSnps)){
    if($$snppipelineSnps{$truePos}){
      $TP++;  # True positive: correct SNP was found
    } else {
      $FN++;  # False negative: a real SNP was not found
    }
  }

  # How many SNPs were found that were not real?
  while(my($pos,$ref)=each(%$snppipelineSnps)){
    if($$realSnps{$pos}){
      # This is a true positive and was already counted in the previous loop.
    } else {
      $FP++;  # A SNP was found but is not in the real set of SNPs.
    }
  }

  return ($TP,$TN,$FP,$FN);
}

sub usage{
  "Compares a Snp-Pipeline run to a simulated dataset
  Usage: $0 --snpPipeline projdirectory --ttr treetoreadsdirectory
  "
}
