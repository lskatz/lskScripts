#!/usr/bin/env perl
# Author: Lee Katz
# Safely remove a taxon from a tree without incurring singleton nodes

use strict;
use warnings;
use Bio::TreeIO;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;

sub logmsg{local $0=basename($0); print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tree=s)) or die $!;
  
  die usage() if($$settings{help});
  $$settings{tree} || die "ERROR: no tree was given:\n".usage();
  
  my @remove=@ARGV;
  die "ERROR: need to remove at least taxon\n".usage() if(@remove < 1);

  my $treeObj=safeRemove($$settings{tree},\@remove,$settings);

  my $out=Bio::TreeIO->new(-format=>"newick");
  $out->write_tree($treeObj);
  print "\n"; # just because newick files don't have newlines for some reason

  return 0;
}

sub safeRemove{
  my($tree,$remove,$settings)=@_;

  my $treeObj=Bio::TreeIO->new(-file=>$tree)->next_tree;

  $treeObj=removeTaxa($treeObj,$remove,$settings);
  $treeObj=removeSingletons($treeObj,$settings);
  
  return $treeObj;
}

sub removeTaxa{
  my($tree,$remove,$settings)=@_;

  my %taxaInTree=();
  my @taxaInTree=$tree->get_leaf_nodes;
  $taxaInTree{$_->id}=1 for(@taxaInTree);

  for my $taxon(@$remove){
    die "ERROR: taxon $taxon does not exist in the tree!" if(!$taxaInTree{$taxon});
    my $bool=$tree->remove_Node($taxon);
  }

  return $tree;
}

# This is actually the contract_linear_paths function
sub removeSingletons{
  my($tree,$settings)=@_;

  $tree->contract_linear_paths(1);

  return $tree;
}

sub usage{
  local $0=basename($0);
  "$0: Removes a taxon from a tree
  Usage: $0 --tree tree.dnd taxon1 [taxon2...]
  "
}
