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
  
  return $treeObj;
}

sub removeTaxa{
  my($tree,$remove,$settings)=@_;

  my %leaf_node_id=();
  my %ancestor_node=();
  my @ancestor_node=();
  my @node = $tree->get_nodes;
  for my $node(@node){
    if($node->is_Leaf()){
      $leaf_node_id{$node->id}=1;
    } else {
      push(@ancestor_node, $node);
      $ancestor_node{$node}=1;
    }
  }

  for my $taxon(@$remove){
    die "ERROR: taxon $taxon does not exist in the tree!" if(!$leaf_node_id{$taxon});
    my $safely_removed=$tree->remove_Node($taxon);
    if(!$safely_removed){
      die "ERROR: could not remove $taxon safely";
    }
  }

  $tree->contract_linear_paths(1);

  # Now remove all nodes that were ancestors but are now leaves
  #my $nodes_were_removed=1;
  #while($nodes_were_removed){
  #  $nodes_were_removed = removeUselessNewLeafNodes($tree,\@ancestor_node);
  #  last;
  #}

  return $tree;
}

sub removeUselessNewLeafNodes{
  my($tree,$ancestor_node)=@_;

  my $nodesRemovedCounter=0;
  for my $node(@$ancestor_node){
    # If an ancestor node is now a leaf, prune it too
    if($node->is_Leaf()){
      my $safely_removed=$tree->remove_Node($node);
      if(!$safely_removed){
        die "ERROR: could not remove a useless ancestor node safely";
      }
      $nodesRemovedCounter++;
    }
  }
  return $nodesRemovedCounter;
}

sub usage{
  local $0=basename($0);
  "$0: Removes a taxon from a tree
  Usage: $0 --tree tree.dnd taxon1 [taxon2...]
  "
}
