#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

use version 0.77;
our $VERSION = '0.1.1';

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help minpercent|min_percent|min-percent=f)) or die $!;
  $$settings{minpercent} ||= 25.0;
  usage() if(!@ARGV || $$settings{help});

  # print header first
  print join("\t", qw(
    file
    BestGuess GuessPercent
    Rank
    Contaminant ContaminantPercent
  )) . "\n";
    
  # Figure out any contamination for each argument
  for my $taxfile(@ARGV){
    my $guesses = guessTaxon($taxfile, $settings);

    print join("\t",
      $taxfile,
      $$guesses{guess}{taxname},
      $$guesses{guess}{percent},
      $$guesses{guess}{rank},
      $$guesses{contaminant}{taxname},
      $$guesses{contaminant}{percent}
    ) . "\n";
  }

  return 0;
}

sub guessTaxon{
  my($taxfile,$settings)=@_;

  # If sn_kraken did not complete, a file will not be present
  if(!-e $taxfile){
    logmsg "WARNING: kraken report not found at $taxfile. It is possible that sn_kraken.pl or kraken was not run on this sample yet.";
    return {}; # return empty hash because that is the var type expected
  }

  my %bestGuess;
  my @header = qw(percent classifiedReads specificClassifiedReads rank taxid taxname);
  open(TAXONOMY, '<', $taxfile) or die "ERROR: could not open $taxfile for reading: $!";
  while(my $taxline = <TAXONOMY>){
    # trim whitespace
    $taxline =~ s/^\s+|\s+$//g;

    # Split fields on tab, but also remove whitespace on fields
    my @F = split(/\s*\t\s*/, $taxline);

    # Name the fields
    my %field;
    @field{@header} = @F;

    # We only care about named ranks like phylum, family, genus, or species
    next if($field{rank} eq '-');

    push(@{ $bestGuess{$field{rank}} }, \%field);

  }
  close TAXONOMY;

  my @sortedRank = qw(S G F O C P K D U);

  # Starting with species, if >min-percent% reads are attributed to
  # a taxon, then guess that taxon.
  my %guessedTaxon;
  RANK:
  for my $rank(reverse @sortedRank){
    $bestGuess{$rank} //= [];
    for my $taxHash(sort {$$b{percent} <=> $$a{percent}} @{ $bestGuess{$rank} }){
      # If we get a good match, starting with highest percentage
      # of reads and going down, then save that into
      # %guessedTaxon and quit the loop.
      if($$taxHash{percent} > $$settings{minpercent}){
        %guessedTaxon = %$taxHash;
      }
    }
  }
  my $rank = $guessedTaxon{rank};

  # Are there any competing taxa at that rank?
  # Initialize the competing taxon to boolean false numbers
  my %majorConflictingTaxon = (rank=>"", percent => 0, taxid=>0, taxname=>".", specificClassifiedReads=>0, classifiedReads=>0);
  $bestGuess{$rank} //= [];
  for my $taxHash(sort {$$b{percent} <=> $$a{percent}} @{ $bestGuess{$rank} }){
    # Skip comparing against self, by comparing taxids
    next if($$taxHash{taxid} == $guessedTaxon{taxid});

    # Accept the conflicting taxon if it's over 1% and if
    # it is the biggest contaminant.
    if($$taxHash{percent} > 1 && $$taxHash{percent} > $majorConflictingTaxon{percent}){
      %majorConflictingTaxon = %$taxHash;
    }
  }

  # Return the best guess and the major conflict
  return {guess=>\%guessedTaxon, contaminant=>\%majorConflictingTaxon};
}


sub usage{
  print "$0: Reads a Kraken report and prints a table of what
  the sample appears to be and then the majority contaminant.

  Usage: $0 [options] kraken.report another.report
  --min-percent What percent of reads have to be attributed
                to a taxon before presuming it as the taxon
                sequenced? Default: 25.0
  --help        This useful help menu
  \n";
  exit 0;
}
