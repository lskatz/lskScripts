#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  die usage() if($$settings{help} || !@ARGV);

  my(@infile) = @ARGV;

  for my $file(@infile){
    logmsg "Reading $file";
    my $alleles = readAlleles($file, $settings);
    my $distances=distances($alleles, $settings);

    my @name = sort(keys(%$distances));
    for(my $i=0;$i<@name;$i++){
      for(my $j=$i+1;$j<@name;$j++){
        my($nameI,$nameJ) = sort ($name[$i],$name[$j]);
        next if($name[$i] !~ /_100/ && $name[$j] !~ /_100/);
        my $distance = $$distances{$name[$i]}{$name[$j]} // "UNKNOWN";
        print join("\t", $name[$i], $name[$j], $distance)."\n";
      }
    }
  }

  return 0;
}

sub readAlleles{
  my($file, $settings) = @_;

  my %allele;

  open(my $fh, $file) or die "ERROR reading $file: $!";
  my $header = <$fh>;
  my @header = split(/\t/, $header);
  shift(@header); # assume the first column is the key and disregard it otherwise
  while(<$fh>){
    chomp;
    my @F = split /\t/;
    my $name = shift(@F);
    my %F;
    @F{@header} = @F;

    $allele{$name} = \%F;
  }
  close $fh;

  return \%allele;
}

sub distances{
  my($alleles, $settings) = @_;

  my %distance;

  my @name = sort keys %$alleles;
  my @allele=sort keys %{$$alleles{$name[0]}};
  my $numNames = @name;
  my $numAlleles=@allele;

  for(my $i=0; $i<$numNames; $i++){
    for(my $j=$i+1; $j<$numNames; $j++){
      # Sort the names to avoid having to make sure that 
      # each "vice versa" distance doesn't have to be
      # calculated.
      my($nameI, $nameJ) = sort($name[$i], $name[$j]);
      # initialize the distances
      $distance{$nameI}{$nameJ} = 0;

      for(my $k=0; $k<$numAlleles; $k++){
        my $alleleI = $$alleles{$nameI}{$allele[$k]};
        my $alleleJ = $$alleles{$nameJ}{$allele[$k]};

        if($alleleI eq '?' || $alleleJ eq '?'){
          next;
        }

        if($alleleI ne $alleleJ){
          $distance{$nameI}{$nameJ}++;
        }
      }
    }
  }

  return \%distance;
}

sub usage{
  "$0: Calculates pairwise distances between genomes in a bionumerics MLST export
  Usage: $0 [options] mlst.tsv > pairwise.tsv
  --help   This useful help menu
  "
}
