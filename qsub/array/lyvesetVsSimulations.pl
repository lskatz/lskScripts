#!/usr/bin/env perl
# Uses a tree to reads directory and a Lyve-SET run to determine
# true positives, false negatives, false positives
#

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help lyveset=s ttr=s));
  die usage() if($$settings{help});
  die "ERROR: need Lyve-SET project name\n".usage() if(!$$settings{lyveset});
  die "ERROR: need TreeToReads project name\n".usage() if(!$$settings{ttr});

  # where are the Lyve-SET SNPs?
  my $lyveSetSnps=lyveSetSnps($$settings{lyveset},$settings);

  # where are the real SNPs?
  my $realSnps=ttrSites($$settings{ttr},$settings);
  
  # compare
  my($TP,$TN,$FP,$FN)=compareSnps($lyveSetSnps,$realSnps,$settings);

  # Report
  print join("\t",qw(L-S TTR TP TN FP FN))."\n";
  print join("\t",$$settings{lyveset},$$settings{ttr},$TP, $TN, $FP, $FN)."\n";
}

sub lyveSetSnps{
  my($proj,$settings)=@_;
  my %pos;

  my $matrix="$proj/msa/out.snpmatrix.tsv";
  open(SNPMATRIX,"<",$matrix) or die "ERROR: could not open $matrix for reading: $!";
  while(<SNPMATRIX>){
    next if(/^#/);
    chomp;
    my($chr,$pos,$ref,@alt)=split /\t/;
    $pos{$pos}=$ref;
  }
  close SNPMATRIX;

  return \%pos;
}

sub ttrSites{
  my($proj,$settings)=@_;
  my %pos;

  my $matrix="$proj/mutsites.txt";
  open(TTRMATRIX,"<",$matrix) or die "ERROR: could not open $matrix for reading: $!";
  while(<TTRMATRIX>){
    chomp;
    my($pos)=split /\s+/;
    $pos{$pos}=1
  }
  close TTRMATRIX;
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
  my ($lyveSetSnps,$realSnps,$settings)=@_;

  # Initialize counts to zero
  my($TP,$TN,$FP,$FN)=split(//,"0" x 4);
  
  # How many of the real SNPs were found?
  while(my($truePos,$trueRef)=each(%$realSnps)){
    if($$lyveSetSnps{$truePos}){
      $TP++;  # True positive: correct SNP was found
    } else {
      $FN++;  # False negative: a real SNP was not found
    }
  }

  # How many SNPs were found that were not real?
  while(my($pos,$ref)=each(%$lyveSetSnps)){
    if($$realSnps{$pos}){
      # This is a true positive and was already counted in the previous loop.
    } else {
      $FP++;  # A SNP was found but is not in the real set of SNPs.
    }
  }

  return ($TP,$TN,$FP,$FN);
}

sub usage{
  "Compares a Lyve-SET run to a simulated dataset
  Usage: $0 --lyveset projdirectory --ttr treetoreadsdirectory
  "
}
