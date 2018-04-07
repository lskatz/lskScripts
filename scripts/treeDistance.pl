#!/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/fileparse basename dirname/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help metric=s lambda=f symmetric)) or die $!;
  die usage() if($$settings{help});
  $$settings{metric}||="kendall";
  $$settings{metric}=lc($$settings{metric});
  $$settings{lambda}||=0.5;

  my @tree=@ARGV;
  my $numTrees=@tree;
  die "ERROR: need at least two trees\n".usage() if($numTrees < 2);
  # Check for file existence
  for(@tree){
    die "ERROR: cannot find $_" if(!-e $_);
  }

  allDistances(\@tree,$settings);

}

# Get all distances
sub allDistances{
  my($tree,$settings)=@_;
  my $numTrees=@$tree;

  # Turn trees to tree objects for speed
  my @treeObj=();
  for(my $i=0;$i<$numTrees;$i++){
    $treeObj[$i]=Bio::TreeIO->new(-file=>$$tree[$i],-format=>"newick")->next_tree;
    reroot($treeObj[$i],$settings);
  }

  # Print header
  print join("\t",qw(FILE1 FILE2 DISTANCE))."\n";

  # Get distances
  for(my $i=0; $i<$numTrees;$i++){
    my $jInit=$i+1;
    if($$settings{symmetric}){
      $jInit=0;
    }
    for(my $j=$jInit;$j<$numTrees;$j++){
      my $distance;
      if($$settings{metric} eq 'kendall'){
        $distance=kendallDistance($treeObj[$i],$treeObj[$j],$settings);
      } else {
        die "ERROR: I do not understand distance metric $$settings{metric}";
      }

      $distance=sprintf("%0.9f",$distance);
      print join("\t",$$tree[$i],$$tree[$j],$distance)."\n";
    }
  }

  return 0;
}

sub kendallDistance{
  my($tree1,$tree2,$settings)=@_;

  # Get all the IDs in the tree. However since both trees
  # are supposed to have the same leaf IDs, only use
  # the array @node1 and not @node2 so that ordering and
  # counts can be correct.
  my @node1=sort {$a->id cmp $b->id } grep {$_->is_Leaf} $tree1->get_nodes;
  my @node2=sort {$a->id cmp $b->id } grep {$_->is_Leaf} $tree2->get_nodes;
  my @id1=map {$_->id} @node1;
  my @id2=map {$_->id} @node1;
  my $root1=$tree1->get_root_node;
  my $root2=$tree2->get_root_node;
  my $numTaxa=@id1;

  # Make sure @node1 eq @node2
  my $errors="";
  for(my $i=0;$i<$numTaxa;$i++){
    if($id1[$i] ne $id2[$i]){
      warn "ERROR: $id1[$i] is not equal to $id2[$i]\n";
      $errors++;
    }
  }
  die "Found $errors errors in the comparison between trees $tree1 and $tree2" if($errors);

  my($m1,$m2,$M1,$M2);
  # m: count of edges in the path
  # For each pair of nodes, figure out the number of steps from their
  # common ancestor to the root.
  my $topologyEuclideanSum=0;
  for(my $i=0;$i<$numTaxa;$i++){
    for(my $j=$i+1;$j<$numTaxa;$j++){
      # Need the number of steps between the lowest common ancestor between
      # the two nodes and the root.
      my $closestAncNode1=$tree1->get_lca([$node1[$i],$node1[$j]]);
      my $closestAncNode2=$tree2->get_lca([$node1[$i],$node1[$j]]);
      my @lineage1=$tree1->get_lineage_nodes($closestAncNode1);
      my @lineage2=$tree2->get_lineage_nodes($closestAncNode2);

      # lineage needs fixing
      @lineage1=grep{defined($_->id) && $_->id ne ''} @lineage1;
      @lineage2=grep{defined($_->id) && $_->id ne ''} @lineage2;
      #next if(!@lineage1 || !@lineage2);

      $topologyEuclideanSum+=(scalar(@lineage1)-scalar(@lineage2))^2;
      #print join("\t",scalar(@lineage1),scalar(@lineage2),$id1[$i],$id1[$j])."\n";
    }
    #print "\n";
  }
  $topologyEuclideanSum=sqrt($topologyEuclideanSum);

  # M: distance of path. In bioperl language, this is the height.
  my $heightEuclideanSum=0;
  for(my $i=0;$i<$numTaxa;$i++){
    for(my $j=$i+1;$j<$numTaxa;$j++){
      my $closestAncNode1=$tree1->get_lca([$node1[$i],$node1[$j]]);
      my $closestAncNode2=$tree2->get_lca([$node1[$i],$node1[$j]]);
      my $height1=$closestAncNode1->height;
      my $height2=$closestAncNode2->height;
      
      $heightEuclideanSum+=($height1-$height2)^2;
    }
  }
  $heightEuclideanSum=sqrt($heightEuclideanSum);

  return abs(($$settings{lambda}*$topologyEuclideanSum) - (1-$$settings{lambda})*$heightEuclideanSum);
}

=cut
    my @node=$treeObj[$i]->get_nodes();
    for(my $k=0;$k<@node;$k++){
      for(my $j=$i+1;$j<@node;$j++){
        print "$k,$j ". $treeObj[$i]->distance([$node[$k],$node[$j]])."\n";
      }
    }
=cut

sub reroot{
  my($tree,$settings)=@_;

  # Set a default outgroup before looking at the rest of the nodes.
  my @node=$tree->get_nodes;
  my $numNodes=@node;
  my $outgroup=$node[0];
  my $longest=$node[0]->branch_length || 0;

  for(my $i=0;$i<$numNodes;$i++){
    for(my $j=$i+1;$j<$numNodes;$j++){
      if($tree->distance([$node[$i],$node[$j]]) > $longest){
        $longest=$tree->distance([$node[$i],$node[$j]]);
        $outgroup=$node[$i];
      }
    }
  }
  $tree->reroot_at_midpoint($outgroup,'MidpointRoot');
  return $outgroup;

  # Need to make sure the tree was rerooted correctly
  my $out=Bio::TreeIO->new(-format=>'tabtree');
  $out->write_tree($tree);
  die;
}

sub usage{
  local $0=basename $0;
  "$0: Finds the distance between two or more trees
  Usage: $0 1.dnd 2.dnd [3.dnd ...]
  --metric     Kendall  Currently, only Kendall is supported
  --symmetric           See if the distance between trees 1
                        and 2 is the same as 2 and 1.
  --lambda     0.5      The kendall coefficient weighing
                        topology vs branch lengths. Values
                        can be zero to one.
  "
}
