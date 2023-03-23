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
  GetOptions($settings,qw(help regex=s essence)) or die $!;
  die usage() if($$settings{help} || !@ARGV);

  my $regex = $$settings{regex} || die "ERROR: --regex is required";

  my ($resultsDir) = @ARGV;
  my $alleles = readAlleles($resultsDir, $settings);
  printDiffs($alleles, $regex, $settings);

  return 0;
}

sub printDiffs{
  my($alleles, $regex, $settings) = @_;

  my @comparedFasta = sort grep {/$regex/} keys(%$alleles);
  my @comparedAllele= sort keys( %{$$alleles{$comparedFasta[0]}} );
  my $numAlleles = @comparedAllele;

  logmsg "Comparing these assemblies: @comparedFasta";

  print join("\t", qw(ref query locus refAllele queryAllele))."\n";

  for(my $i=1; $i<@comparedFasta; $i++){
    for(my $j=0; $j<$numAlleles; $j++){
      my $refAllele   = $$alleles{$comparedFasta[0]}{$comparedAllele[$j]};
      my $queryAllele = $$alleles{$comparedFasta[$i]}{$comparedAllele[$j]};
      my $refLocus  = $comparedAllele[$j];

      # Error checking
      $refAllele   //= die "ERROR no allele found for $comparedFasta[0] at $refLocus";
      $queryAllele //= die "ERROR no allele found for $comparedFasta[$i] at $refLocus";

      if($$settings{essence}){
        # Skip lines where the locus was not found in either query or ref
        next if("$queryAllele$refAllele" =~ /LNF/);
        next if("$queryAllele$refAllele" =~ /PLOT/);
        next if("$queryAllele$refAllele" =~ /NIPH|NIPHEM/);
        next if("$queryAllele$refAllele" =~ /LOTSC/);
        next if("$queryAllele$refAllele" =~ /A[LS]M/);

        # Remove *
        $refAllele   =~ s/\*//g;
        $queryAllele =~ s/\*//g;

        # Remove INF
        $refAllele   =~ s/INF-//;
        $queryAllele =~ s/INF-//;
      }

      if($refAllele ne $queryAllele){
        print join("\t", $comparedFasta[0], $comparedFasta[$i], $refLocus, $refAllele, $queryAllele);
        print "\n";
      }
    }
  }
}

sub readAlleles{
  my($resultsDir, $settings) = @_;

  my %alleles;

  my $alleles = "$resultsDir/results_alleles.tsv";
  open(my $fh, $alleles) or die "ERROR: could not read $alleles: $!";

  # So tired of using multiple lines to get the header.
  # Why not get it in one go.
  chomp(my @header = split(/\t/, scalar(<$fh>)));
  # remove the header FILE
  my $FILE = shift(@header);
  die "INTERNAL ERROR: \$FILE is not literally FILE" if($FILE ne "FILE");

  while(<$fh>){
    chomp;
    my ($filename, @allele) = split(/\t/, $_);
    
    my %F;
    @F{@header} = @allele;

    $alleles{$filename} = \%F;
  }
  close $fh;

  return \%alleles;
}


sub usage{
  "$0: Find allelic differences from chewie results dir
  Usage: $0 [options] resultsdir
  --regex  Compare assemblies in the results with this regex
           Example to compare Campy_shovill_SRR13093829_skesa.fasta and Campy_shovill_SRR13093829_spades.fasta:
           --regex 'Campy_shovill_SRR13093829_(skesa|spades)'
  --essence  
           Remove * and INF- from allele names before comparing
           Remove lines with LNF, PLOT, NIPH, NIPHEM
           => get just actual differences
  --help   This useful help menu
  "
}
