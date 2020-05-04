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

  my @file = @ARGV;
  my $numFiles = @file;
  my %dist;
  for my $file(@file){
    $dist{$file} = readPhylip($file, $settings);
  }

  my @taxon = keys(%{ $dist{$file[0]} });
  my $numTaxa = @taxon;

  print join("\t", "taxon1", "taxon2", @file)."\n";
  for(my $i=0;$i<$numTaxa;$i++){
    for(my $j=0;$j<$numTaxa;$j++){
      print $taxon[$i]."\t".$taxon[$j];
      for(my $k=0;$k<$numFiles;$k++){
        my $singleDist = $dist{$file[$k]}{$taxon[$i]}{$taxon[$j]};
        if(!defined($singleDist)){
          $singleDist = $dist{$file[$k]}{$taxon[$j]}{$taxon[$i]};
        }
        if(!defined($singleDist)){
          die "ERROR: distance not found:\n".Dumper [$k,$file[$k]], [$i,$taxon[$i]], [$j, $taxon[$j]];
        }
        print "\t".$singleDist;
      }
      print "\n";
    }
  }

  return 0;
}

sub readPhylip{
  my($file, $settings)=@_;
  
  my %dist;
  my %distArr;
  my @taxon;
  open(my $fh, $file) or die "ERROR: could not read $file: $!";
  my $numTaxa = <$fh>;
  $numTaxa =~ s/^\s+|\s+$//g;
  while(<$fh>){
    chomp;
    my($taxon, @dist) = split(/\s+/, $_);
    $distArr{$taxon} = \@dist;
    push(@taxon, $taxon);
  }

  my $actualNumTaxa = @taxon;
  if($actualNumTaxa != $numTaxa){
    die "ERROR: in $file, reported number of taxa does not match number of taxa found";
  }
  
  # Now that we know all the taxa in the list, go back
  # and fill in the 2d hash
  while(my($refTaxon, $distances) = each(%distArr)){
    for(my $i=0;$i<$numTaxa;$i++){
      my $queryTaxon = $taxon[$i];
      $dist{$refTaxon}{$queryTaxon} = $$distances[$i];
    }
  }

  return \%dist;
}

sub usage{
  "
  $0: changes phylip distance files to a single tall/skinny format
  Usage: $0 [options] file1.phylip [file2.phylip...]
  --help   This useful help menu

  "
}
