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
  GetOptions($settings,qw(help root=s));
  $$settings{root}||="";
  die usage() if($$settings{help});

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
    reroot($tree,$settings) if($$settings{root});
    for my $node($tree->get_nodes(-order=>"depth")){ # other choice: "breadth"
      next if(!$node->is_Leaf);
      my $id=$node->id;
      $id=~s/^'|'$//g; # remove single quotes at beginning/end that bioperl adds
      print "$id\n";
    }
  }
  $in->close;
}

sub reroot{
  my($tree,$settings)=@_;
  # TODO validate that {root} eq 'midpoint' or a node name
  if($$settings{root} =~/^longest$/i){
    # converge on a longest branch
    _rerootLongestBranch($tree,$settings) for(1..10);
  } else {
    # Reroot on whichever node matches.
    # NOTE: the name is not validated here and so 
    # rerooting might not happen if the ID is not found.
    for my $node($tree->get_leaf_nodes){
      $tree->reroot($node) if($node->id eq $$settings{root});
    }
  }
}

sub _rerootLongestBranch{
  my($tree,$settings)=@_;

  my @node=$tree->get_nodes;
  my $outgroup=$node[0];
  my $longest=$node[0]->branch_length || 0;

  for(my $i=1;$i<@node;$i++){
    if($node[$i]->branch_length() > $longest){
      $longest=$node[$i]->branch_length;
      $outgroup=$node[$i];
      #print join("\t",$outgroup->id,$longest)."\n";
    }
  }

  $tree->reroot($outgroup);
}

sub usage{
  "$0: determine the phylogenetic order from a tree file
  Usage: $0 tree.dnd
  --root longest  Reroot the tree at the longest branch.  If you supply a taxon ID then it will root on the branch leading to it instead.
  "
}
