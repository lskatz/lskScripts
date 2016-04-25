#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;

  my @tree=@ARGV;

  print join("\t",qw(File numLeaves numNodes))."\n";
  for my $t(@tree){
    treeInfo($t,$settings);
  }
  
  return 0;
}

sub treeInfo{
  my($tree,$settings)=@_;

  my $t=Bio::TreeIO->new(-file=>$tree)->next_tree;
  my @leaf=map{$_->id} $t->get_leaf_nodes;
  my @node=map{$_->id} $t->get_nodes;

  print join("\t",$tree,scalar(@leaf),scalar(@node))."\n";
}

sub usage{
  $0=basename $0;
  "$0: gives information about a phylogeny
  Usage: $0 tree.dnd [tree2.dnd...]
  "
}
