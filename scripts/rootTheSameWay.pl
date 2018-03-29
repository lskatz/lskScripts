#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;
use List::Util qw/sum/;
use Statistics::Descriptive;
use File::Temp qw/tempdir tempfile/;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;

  die usage() if(!@ARGV || $$settings{help});

  my($ref,@query)=@ARGV;
  my $treeout=Bio::TreeIO->new(-format=>"newick");

  my $refObject        = Bio::TreeIO->new(-file=>$ref)->next_tree;
  my $refTreeOutgroup  = treeOutgroup($refObject,$settings);
  my @outgroupIds      = map{$_->id} grep{$_->is_Leaf} $refTreeOutgroup->get_all_Descendents($refTreeOutgroup);
  my @leaves           = map{$_->id} grep{$_->is_Leaf} $refObject->get_nodes;

  # For each query, reroot and print
  for my $q(@query){
    my $queryObject      = Bio::TreeIO->new(-file=>$q)->next_tree;
    my $queryTreeOutgroup= findNode([@outgroupIds,@leaves],$queryObject,$settings);
    logmsg "Found the correct node to root the tree";
    $queryObject->reroot($queryTreeOutgroup);

    my @sortedQueryId = sort {$a cmp $b} map {$_->id} grep {$_->is_Leaf} 
      ($queryObject->get_root_node->each_Descendent)[0]->get_all_Descendents;
    logmsg "$sortedQueryId[0]...",scalar(@sortedQueryId),scalar(@outgroupIds);

    $treeout->write_tree($queryObject);
  }

  return 0;
}

# Given the way the tree is rooted, return which node
# is the outgroup, defined by the fewest number of leaves
# of the first descendent of the root.
sub treeOutgroup{
  my($treeObj,$settings)=@_;

  # Get all direct descendents of the root node,
  # sorted by number of descendents.
  my @secondaryNode=sort{
    my @taxaA=sort {$a->id cmp $b->id} grep{$_->is_Leaf} $a->get_all_Descendents;
    my @taxaB=sort {$a->id cmp $b->id} grep{$_->is_Leaf} $b->get_all_Descendents;
    my $numA=scalar(@taxaA);
    my $numB=scalar(@taxaB);

    return 0 if($numA==0 && $numB==0);

    if($numA != $numB){
      return $numA <=> $numB;
    }

    # If they have the same number of taxa, go with
    # alphabetical order to keep a stable sort
    for(my $i=0;$i<$numA;$i++){
      if($taxaA[$i] ne $taxaB[$i]){
        return $taxaA[$i] cmp $taxaB[$i];
      }
    }
    
    die "INTERNAL ERROR: somehow this tree bifurcates into clades with equally named leaves???";
  } $treeObj->get_root_node->each_Descendent;
  for my $node(@secondaryNode){
    my $numLeaves=scalar(grep{$_->is_Leaf} $node->get_all_Descendents);
    if($numLeaves > 0){
      return $node;
    }
  }
  
  return $secondaryNode[0];
}

sub findNode{
  my($leafIds,$tree,$settings)=@_;

  my @sortedLeafIds = sort{$a cmp $b} @$leafIds;
  my $numLeaves=@$leafIds;

  # Look for the node first by starting at one leaf and 
  # moving up the tree by ancestor.
  for(my $nodeIndex=0;$nodeIndex<$numLeaves;$nodeIndex++){
    my @starterNode=$tree->find_node(-id=>$$leafIds[$nodeIndex]);
    for(my $i=0;$i<@starterNode;$i++){
      my $node=$starterNode[$i];
      if($node->each_Descendent < 2){
        next;
      }
      # Go up the chain until we have at least the right number
      # of decendent nodes.
      my @sortedQueryId;
      do{
        @sortedQueryId = sort {$a cmp $b} map {$_->id} grep {$_->is_Leaf} $node->get_all_Descendents;
        if("@sortedLeafIds" eq "@sortedQueryId"){
          $node->id("theRoot");
          return $node;
        }
        $node=$node->ancestor;
      } while(@sortedQueryId < $numLeaves);

      #while(scalar($node->get_all_Descendents) < $numLeaves){
      #  $node=$node->ancestor;
      #  my @sortedQueryId = sort {$a cmp $b} map {$_->id} grep {$_->is_Leaf} $node->get_all_Descendents;
      #  if("@sortedLeafIds" eq "@sortedQueryId"){
      #    return $node;
      #  }
      #}
    }
  }
  die "Could not find common root with the heuristic.";


  # If all else fails, just try every node and reroot

  # look at all nodes and find the one with the correct descendents
  #for my $node($tree->get_nodes)
}


sub usage{
  "$0: Roots a query tree according to a reference tree
  Usage: $0 ref.dnd query.dnd > query.rerooted.dnd
  "
}

