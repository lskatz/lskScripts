#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Perl;
use Bio::TreeIO;
use Getopt::Long;
use Data::Dumper;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help));
  exit usage() if($$settings{help});

  my(@tree)=@ARGV;
  die "ERROR: need a tree file!\n".usage() if(!@tree);

  for my $tree(@tree){
    printPhylogeneticOrder($tree,$settings);
  }
  return 0;
}

sub printPhylogeneticOrder{
  my($tree,$settings)=@_;

  my $numtaxa=0;
  my $in=Bio::TreeIO->new(-file=>$tree);
  while(my $tree=$in->next_tree){
    #reroot($tree,$settings);
    for my $node($tree->get_nodes(-order=>"depth")){ # other choice: "breadth"
      next if(!$node->is_Leaf);
      print $node->id."\n";
    }
  }
  $in->close;
}

# TODO use $node->branch_length to determine the longest branch length for rooting
sub reroot{
  my($tree,$settings)=@_;
  my @node=$tree->get_leaf_nodes;
  my %distance;
  my $numnodes=@node;
  my $largestDistance=0;
  for(my $i=0;$i<$numnodes;$i++){
    my $id1=$node[$i]->id;
    for(my $j=0;$j<$numnodes;$j++){
      my $id2=$node[$j]->id;
      $distance{$id1}{$id2}=$tree->distance(-nodes=>[$node[$i],$node[$j]])
    }
  }
  die Dumper \%distance;
  #tree->reroot($farthestNode);
}

sub usage{
  "$0: determine the phylogenetic order from a tree file
  Usage: $0 tree.dnd
  "
}
