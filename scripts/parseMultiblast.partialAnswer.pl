#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;

my %query=();
my %hit=();
my $section="";
while(<>){ # read line by line
  if(/Table of genes/){
    $section="Table of genes";
    while(<>){
      # The section ends on a blank line
      if(/^\s*$/){
        last;
      }

      chomp; # remove whitespace
      my($locus, $start, $stop, $strand, $annotation)=split(/\s+/);
      $query{$locus}={
        start     =>$start,
        stop      =>$stop,
        strand    =>$strand,
        annotation=>$annotation,
      };
      
    }
  }

  # Redundant with the details section
  elsif(/Significant hits/){
    while(<>){
      if(/Details/){
        $section="Details";
        parseDetailsSection(\%hit, \%query);
        last;
      }
    }
  }
}

print Dumper \%query;

sub parseDetailsSection{
  my($hit,$query)=@_;

  my $currentHit="";
  while(<>){
    if(/^\s*$/){
      next;
    }

    chomp;
    if(/\d+\.\s+(\S+)/){
      $currentHit=$1;
      # Source is on the next line.
      my $source = scalar(<>);
      chomp $source;
      $source=~s/^Source: //;
      $$hit{$currentHit}{source}=$source;
    }

    # A regex double meaning with elipses: a commentary
    # that this is too wordy but it also simply matches 
    # on at least three characters.
    # Also: this perl comment is wordy.
    if(/^Number of proteins with BLAST hits...+(\d+)/){
      $$hit{$currentHit}{numhits}=$1;
      $$hit{$currentHit}{multiblastscore}=scalar(<>);
      $$hit{$currentHit}{multiblastscore}=~s/\s+|\D+//g; # trim and remove non-digits
      $$hit{$currentHit}{blastscore}=scalar(<>)i
      $$hit{$currentHit}{blastscore}=~s/\s+|\D+//g; # trim and remove non-digits
    }
    elsif(/Table of genes.../){
      while(<>){
        if(/^\s*$/){
          last;
        }
        chomp;
        my($locus, $start, $stop, $strand, $annotation)=split(/\s+/);
        $$hit{$currentHit}{query}{$locus}={
          start     => $start,
          stop      => $stop,
          strand    => $strand,
          annotation=> $annotation,
        };
    }

  }
}
