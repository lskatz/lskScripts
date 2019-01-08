#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Bio::TreeIO;
use Bio::Tree::Statistics;

my $guideTree; # will be the first tree
my @bs_tree;
for my $file(@ARGV){
  my $in=Bio::TreeIO->new(-file=>$file);
  while(my $tree=$in->next_tree){
    if(!$guideTree){
      $guideTree = $tree;
      next;
    }
    push(@bs_tree, $tree);
  }
}

my $biostat = Bio::Tree::Statistics->new();
my $bsTree=$biostat->assess_bootstrap(\@bs_tree,$guideTree);
for my $node($bsTree->get_nodes){
  if(!$node->id){
    $node->id($node->bootstrap);
  }
}
print $bsTree->as_text("newick")."\n";
