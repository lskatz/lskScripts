#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
# Print to screen whatever warning/error in the taxonomy
sub logInvalid{
  my($msg, $settings) = @_;
  $msg = "WARNING: $msg";
  logmsg $msg;
  return 1;
}

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  die usage() if($$settings{help} || @ARGV < 2);

  my($nodes, $names) = @ARGV;

  my $error_counter = validateTaxonomy($nodes, $names, $settings);

  logmsg "$error_counter errors found in $nodes and $names";

  return 0;
}

sub validateTaxonomy{
  my($nodes, $names, $settings) = @_;

  # Read the files into hashes
  logmsg "Reading $names";
  my $taxonNames = readNames($names, $settings);
  logmsg "Reading $nodes";
  my $lineage    = readNodes($nodes, $settings);

  # Start looping through the taxa. Rely on the taxids
  # from the lineages.
  my $error_counter = 0;
  my @taxid = sort {$a<=>$b} keys(%$lineage);
  for my $id(@taxid){

    # Check to make sure that it has a names entry
    my $namesEntry = $$taxonNames{$id};
    if(! $namesEntry){
      $error_counter += logInvalid "no entry for $id";
    }

    # Check to make sure that the names entry has a name class of scientific name
    my $scientificName = $$namesEntry{'scientific name'};
    if(! $scientificName){
      $error_counter += logInvalid "no scientific name found for $id";
    }

    # Check to make sure that there is a parent taxid
    if(! $$lineage{$id}{parent_id}){
      $error_counter += logInvalid "no parent id for $id ($scientificName)";
    }

    # Check to make sure the parent taxid exists in lineage
    my $parent_id = $$lineage{$id}{parent_id};
    if(! $$lineage{$parent_id}){
      $error_counter += logInvalid "Parent ID $parent_id (parent of $id) not found";
    }
  }

  return $error_counter;
}

sub readNames{
  my($names, $settings) = @_;

  my %name;

  # For file IO, change the character to newline to "\t|\n"
  # to match the file specifications.
  local $/ = "\t|\n";
  open(my $fh, $names) or die "ERROR: could not read names file $names: $!";
  my @header = qw(tax_id name_txt unique_name name_class);
  while(<$fh>){
    chomp;
    my @F = split(/\t\|\t/, $_);
    my %F;
    @F{@header} = @F;

    $name{$F{tax_id}}{$F{name_class}} = \%F;
  }
  close $fh;

  return \%name;
}

sub readNodes{
  my($nodes, $settings) = @_;

  my %node;

  # For file IO, change the character to newline to "\t|\n"
  # to match the file specifications.
  local $/ = "\t|\n";
  open(my $fh, $nodes) or die "ERROR: could not read nodes file $nodes: $!";
  my @header = qw(tax_id parent_id rank embl division_id
    inherited_div_flag genetic_code_id inherited_MCG_flag
    genbank_hidden_flag hidden_subtree_root_flag comments
  );
  while(<$fh>){
    chomp;
    my @F = split(/\t\|\t/, $_);
    my %F;
    @F{@header} = @F;

    $node{$F{tax_id}} = \%F;
  }
  close $fh;

  return \%node;
}

sub usage{
  "$0: validates a taxononomic scheme with nodes.dmp and names.dmp
  Usage: $0 [options] nodes.dmp names.dmp
  --help   This useful help menu
  "
}
