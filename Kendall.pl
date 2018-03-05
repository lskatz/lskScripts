#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;
use List::Util qw/sum/;

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help lambda=f alreadyrooted|rooted)) or die $!;

  $$settings{lambda}||=0;

  die usage() if(!@ARGV || $$settings{help});

  my($tree1,$tree2)=@ARGV;

  validateTrees($tree1,$tree2,$settings);

  my $value=kc($tree1,$tree2,$$settings{lambda},$settings);
  print $value."\n";

  return 0;
}

# Find the Kendall-Colijn value of a tree comparison
sub kc{
  my($tree1,$tree2,$lambda,$settings)=@_;

  
  my ($dist1,$numBranches1) = distanceMatrix($tree1,$settings);
  my ($dist2,$numBranches2) = distanceMatrix($tree2,$settings);

  my @leafName = map{$_->id} grep {$_->is_Leaf} Bio::TreeIO->new(-file=>$tree1)->next_tree->get_nodes;
  my $numLeaves=@leafName;

  my @topologySquares;
  my @lengthsSquares;
  for(my $i=0;$i<$numLeaves;$i++){
    # Take care of distance from leaf to ancestor node
    push(@topologySquares,
      ($$numBranches1{$leafName[$i]}{$leafName[$i]} - $$numBranches2{$leafName[$i]}{$leafName[$i]}) ** 2
    );
    push(@lengthsSquares,
      ($$dist1{$leafName[$i]}{$leafName[$i]} - $$dist2{$leafName[$i]}{$leafName[$i]}) ** 2
    );
    for(my $j=$i+1;$j<$numLeaves;$j++){
      push(@topologySquares,
        ($$numBranches1{$leafName[$i]}{$leafName[$j]} - $$numBranches2{$leafName[$i]}{$leafName[$j]}) ** 2
      );
      push(@lengthsSquares,
        ($$dist1{$leafName[$i]}{$leafName[$j]} - $$dist2{$leafName[$i]}{$leafName[$j]}) ** 2
      );
    }
  }

  my $lengthsDistance = sqrt(sum(@lengthsSquares));
  my $topologyDistance= sqrt(sum(@topologySquares));

  my $kc = (1-$lambda) * $topologyDistance + $lambda * $lengthsDistance;
  return $kc;
}

# Find the distance in both branch counts and branch lengths
sub distanceMatrix{
  my($tree,$settings)=@_;
  my $treeObj = Bio::TreeIO->new(-file=>$tree)->next_tree;

  my %distance;
  my %numBranches;

  my @node=grep {$_->is_Leaf} $treeObj->get_nodes;
  my $numNodes=@node;
  for(my $i=0;$i<$numNodes;$i++){
    for(my $j=$i+1;$j<$numNodes;$j++){
      my($dist,$numBranches)=distanceToRoot($treeObj,[$node[$i],$node[$j]]);
      $distance{$node[$i]->id}{$node[$j]->id}=$dist;
      $distance{$node[$j]->id}{$node[$i]->id}=$dist;

      $numBranches{$node[$i]->id}{$node[$j]->id}=$numBranches;
      $numBranches{$node[$j]->id}{$node[$i]->id}=$numBranches;
    }
  }

  # KC also requires the distance from leaf to nearest ancestor node.
  # We can store that as "self vs self"
  for(my $i=0;$i<$numNodes;$i++){
    $distance{$node[$i]->id}{$node[$i]->id}=$node[$i]->branch_length;
    #$numBranches{$node[$i]->id}{$node[$i]->id}=scalar($treeObj->get_lineage_nodes($node[$i]));
    $numBranches{$node[$i]->id}{$node[$i]->id}=1;
  }

  # debug
  #for(my $i=0;$i<$numNodes;$i++){
  #    print join("\t",$node[$i]->id,values(%{$numBranches{$node[$i]->id}}))."\n";
  #}

  return (\%distance,\%numBranches);
}

# Find the distance from an LCA to the root
sub distanceToRoot {
    my ($tree,$nodes) = @_;

    my $num_branches = 0;
    my $branch_length= 0;
    my $lca = $tree->get_lca(@{$nodes});

    unless($lca) { 
        $tree->warn("could not find the lca of supplied nodes; can't find distance either");
        return;
    }

    my $root_node=$tree->get_root_node;
    my $curr_node = $lca;
    while($curr_node ne $root_node){
      $num_branches++;
      $branch_length+=$curr_node->branch_length;
      $curr_node=$curr_node->ancestor;
    }
    return ($branch_length,$num_branches);
}


# For KC, the trees have to be rooted the same way, so test
# the outgroup to make sure it is the same.
# Also test to make sure all leaves have the same name.
sub validateTrees{
  my($tree1,$tree2,$settings)=@_;

  # Test for same root by testing the outgroup leaf names
  if(!$$settings{alreadyrooted}){
    testForOutgroup($tree1,$tree2,$settings);
  }

  testForSameLeaves($tree1,$tree2,$settings);

  return 1;
}

# Test for same leaves
sub testForSameLeaves{
  my($tree1,$tree2,$settings)=@_;
  my $treeObj1 = Bio::TreeIO->new(-file=>$tree1)->next_tree;
  my $treeObj2 = Bio::TreeIO->new(-file=>$tree2)->next_tree;
  my $nodes1_str=join(" ",
    sort{$a cmp $b} map{$_->id} grep{$_->is_Leaf} $treeObj1->get_nodes
  );
  my $nodes2_str=join(" ",
    sort{$a cmp $b} map{$_->id} grep{$_->is_Leaf} $treeObj2->get_nodes
  );
  if($nodes1_str ne $nodes2_str){
    die "ERROR: leaves are not the same in both trees";
  }
  return 1;
}

# Die if the outgroups are not the same.
sub testForOutgroup{
  my($tree1,$tree2,$settings)=@_;
  my $treeObj1 = Bio::TreeIO->new(-file=>$tree1)->next_tree;
  my $treeObj2 = Bio::TreeIO->new(-file=>$tree2)->next_tree;
  my $outgroup1Node = treeOutgroup($treeObj1,$settings);
  my $outgroup2Node = treeOutgroup($treeObj2,$settings);
  my @outgroup1=sort{$a cmp $b} map{$_->id} grep{$_->is_Leaf} $outgroup1Node->get_all_Descendents;
  my @outgroup2=sort{$a cmp $b} map{$_->id} grep{$_->is_Leaf} $outgroup2Node->get_all_Descendents;

  my $outgroup1_str=join(" ",@outgroup1);
  my $outgroup2_str=join(" ",@outgroup2);

  if($outgroup1_str ne $outgroup2_str){
    die "ERROR: outgroups are not the same on the rooted trees:\n"
       .$outgroup1_str."\n"
       .$outgroup2_str."\n";
  }
  return 1;
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


sub usage{
  $0=basename $0;
  "$0: runs the Kendall-Colijn metric
  Usage: $0 tree.dnd tree2.dnd

  --lambda         0  A lambda value to use in the metric.
                      0 gives weight to a topology metric.
                      1 gives weight to branch lengths.
                      Must be between 0 and 1.
  --alreadyrooted     The tree is already rooted; don't try
                      to validate it.
  "
}

