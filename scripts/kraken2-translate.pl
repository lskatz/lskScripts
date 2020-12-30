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
  usage() if($$settings{help} || !@ARGV);

  my($infile) = @ARGV;

  translateKraken2($infile);

  return 0;
}

# Emulate the kraken-translate format with 2-space-delimited columns:
#   count
#   parent node 1
#   parent node 2
#   ...
#   genus
#   genus/species
#   ... further divisions
sub translateKraken2{
  my($infile) = @_;
  
  # An array of parent nodes excluding the current node,
  # going from root to kingdom to phylum all the way down,
  # where the 0th element is almost always going to be root
  my @parent = ();
  my $parentLevel  = 0; # always equivalent to number of spaces x 2

  open(my $fh, "<", $infile) or die "ERROR: could not open $infile: $!";
  while(my $line = <$fh>){
    # left and right whitespace trim
    $line =~ s/^\s+|\s+$//g;

    # fields defined in the manual at
    #   https://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports
    #   * Percentage of reads covered by the clade rooted at this taxon
    #   * Number of reads covered by the clade rooted at this taxon
    #   * Number of reads assigned directly to this taxon
    #   * A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
    #   * NCBI taxonomy ID
    #   * indented scientific name
    my @F = split(/\t/, $line);
    my($percent, $readsUmbrella, $readsSpecific, $rank, $taxid, $nameWithIndentation) = @F;

    # Figure out the level from the number of spaces in the name.
    # Each two spaces is one more level down.
    # There are ways to optimize this step in perl but I just wanted
    # to make it clear what was happening.
    my $childLevel = 0;
    if($nameWithIndentation =~ /^( +)/){
      my $prefixWhitespace = $1;
      while($prefixWhitespace =~ /  /g){
        $childLevel++;
      }
    }

    # Remove the indentation into $name
    my $name = $nameWithIndentation;
    $name =~ s/^\s+//;
    $name =~ s/\s+/_/g; # also remove internal whitespace

    # Cut down the parent nodes to this child level
    # while getting the current node tacked on.
    $parent[$childLevel] = $name;
    
    # Set up the taxa fields with the parent node(s) and the current node.
    my @taxaField = @parent[0..$childLevel];

    print join("\t", $readsSpecific, @taxaField)."\n";
  }
}

sub usage{
  print "$0: changes kraken raw output to a format for ktImportText in Krona
  Usage: $0 [options] kraken.out
  --help   This useful help menu

  Output is tab delimited:
  * count of reads
  * parent node 1
  * ...
  * last child node (usually genus/species)
";
  exit 0;
}
