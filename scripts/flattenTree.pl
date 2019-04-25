#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;
use Scalar::Util qw/looks_like_number/;

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help confidence|bootstrap|min-confidence=f)) or die $!;
  $$settings{confidence}||=0;

  for my $file(@ARGV){
    my $in = Bio::TreeIO->new(-file=>$file);
    while(my $tree = $in->next_tree){
      flattenTree($tree, $$settings{confidence}, $settings);
    }
  }

  return 0;
}

sub flattenTree{
  my($tree, $minConfidence, $settings)=@_;
  
  for my $leaf($tree->get_nodes()){
    next if(!$leaf->is_Leaf());

    my @lineage = ($leaf,reverse($tree->get_lineage_nodes($leaf)));

    for(my $i=0;$i<@lineage;$i++){
      next if(!defined($lineage[$i]->ancestor));

      my $confidence = $lineage[$i]->ancestor->id;
      $confidence //= 0;
      if(!looks_like_number($confidence)){
        next;
      }

      if($confidence < $minConfidence){
        $lineage[$i]->ancestor(
          $lineage[$i]->ancestor->ancestor
        );
      }
    }
  }

  my $numRemoved = 1;
  while($numRemoved > 0){
    $numRemoved = 0;
    # Remove singleton paths
    $tree->contract_linear_paths;
    # Remove dead ancestor nodes
    for my $leaf($tree->get_nodes){
      next if(!$leaf->is_Leaf);

      if(looks_like_number($leaf->id)){
        $tree->remove_Node($leaf);
        $numRemoved++;
      }
    }
  }

  print $tree->as_text('newick')."\n";
}

sub usage{
  $0=basename $0;
  "$0: flattens a tree using node confidence scores
  NOTE: leaves with number-only identifiers will be removed.

  Usage: $0 tree.dnd [tree2.dnd...] > out.dnd
  --confidence   0  Minimum confidence for flattening a tree
  "
}
