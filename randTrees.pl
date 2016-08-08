#!/usr/bin/env perl
use warnings;
use strict;
use Bio::Tree::RandomFactory;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;
use Data::Dumper qw/Dumper/;
use List::Util qw/min max shuffle/;
use List::MoreUtils qw/uniq/;

use threads;
use Thread::Queue;
use threads::shared;

my $numTrees :shared = 0;   # number of trees that have been printed

local $0=basename $0;
sub logmsg { print STDERR "$0: @_\n"}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help numTrees=i numcpus=i force-binary)) or die $!;
  $$settings{numTrees}||=1;
  $$settings{numcpus}||=1;
  $$settings{'force-binary'}||=0;

  my @tree=@ARGV;

  die usage() if($$settings{help} || !@tree);
  
  # read the trees for genome names and branch lengths
  my @branchLength;
  my @sampleName;
  for my $t(@tree){
    my($b,$s)=readSamplesAndBranchs($t,$settings);
    push(@branchLength,@$b);
    push(@sampleName,@$s);
  }
  @sampleName=uniq(@sampleName);
  @branchLength=shuffle(@branchLength);
  my $numBranchLength=@branchLength;

  my @thr;
  for (0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&randTreeWorker,[@sampleName],min(@branchLength),max(@branchLength),$settings);
  }

  # Join all threads safely
  for(@thr){
    $_->join;
  }

  return 0;
}

sub randTreeWorker{
  my($sampleName,$minBranchLength,$maxBranchLength,$settings)=@_;

  # Create as many random trees as requested
  my @printBuffer;
  my $maxTrees=$$settings{numTrees}/$$settings{numcpus} + 1;
  my $factory=Bio::Tree::RandomFactory->new(-taxa=>$sampleName,-maxcount=>$maxTrees);
  while(my $tree=$factory->next_tree){
    if($$settings{'force-binary'}){
      logmsg "not binary: ".$tree->as_text("newick") if(!$tree->is_binary());
      $tree->contract_linear_paths();
      $tree->force_binary();
      die "not binary: ".$tree->as_text("newick") if(!$tree->is_binary());
    }
    # Alter the nodes to my liking
    my $longestBranchLength=0;
    my $longestNode;
    for my $node($tree->get_nodes){
      # Get random branch lengths from the original trees
      my $newBranchLength=$minBranchLength+rand($maxBranchLength-$minBranchLength);
      $node->branch_length($newBranchLength);
      # Give it a high bootstrap because we're not randomizing that
      if(!$node->is_Leaf){
        $node->bootstrap(100);
        $node->id(100);
      }

      if($newBranchLength > $longestBranchLength && $tree->get_root_node ne $node){
        $longestBranchLength=$newBranchLength;
        $longestNode=$node;
      }
    }
    $tree->reroot_at_midpoint($longestNode);
    push(@printBuffer,$tree->as_text("newick")."\n");
  }
  
  # Print
  {
    lock($numTrees);
    for(@printBuffer){
      $numTrees++;
      last if($numTrees > $$settings{numTrees});
      print $_;
    }
  }

  return 0;
}

sub readSamplesAndBranchs{
  my($t,$settings)=@_;
  my(@branchLength,@sampleName);
  my $treein=Bio::TreeIO->new(-file=>$t);
  while(my $tree=$treein->next_tree){
    # Sample names
    my @samples=$tree->get_leaf_nodes;
    for(@samples){
      push(@sampleName,$_->id);
    }
    
    # Branch lengths
    my @node=$tree->get_nodes;
    for my $n(@node){
      next if(!defined($n->branch_length));
      push(@branchLength,$n->branch_length);
    }
  }
  
  return (\@branchLength,\@sampleName);
}

sub usage{
  "$0: uses real trees to create random trees. 
  Min/max of all tree branch lengths will be used.
  A list of all unique taxon names will be used.
  Usage: $0 realtree.dnd [realtree2.dnd...]
  --numTrees           1   How many trees to produce
  --numcpus            1   Cpus to use
  --force-binary           Forces a binary tree: contracts
                           single-descendent nodes;
                           splits multifurcating nodes
                           such that there are exactly two
                           descendents per ancestor node.
  "
}
