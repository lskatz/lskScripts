#!/usr/bin/env perl
# Author: Lee Katz
# Safely remove a taxon from a tree without incurring singleton nodes

use strict;
use warnings;
use Bio::TreeIO;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;

sub logmsg{print STDERR "@_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tree=s)) or die $!;
  
  die usage() if($$settings{help});
  $$settings{tree} || die "ERROR: no tree was given:\n".usage();
  
  my @remove=@ARGV;
  die "ERROR: need to remove at least taxon\n".usage() if(@remove < 1);

  safeRemove($$settings{tree},\@remove,$settings);

  return 0;
}

sub safeRemove{
  my($tree,$remove,$settings)=@_;

  my $treeObj=Bio::TreeIO->new(-file=>$tree)->next_tree;

  removeTaxa($treeObj,$remove,$settings);
  removeSingletons($treeObj,$remove,$settings);

  die;
}

sub removeTaxa{
  my($tree,$remove,$settings)=@_;

  my %taxaInTree=();
  my @taxaInTree=$tree->get_leaf_nodes;
  $taxaInTree{$_->id}=1 for(@taxaInTree);

  for my $taxon(@$remove){
    die "ERROR: taxon $taxon does not exist in the tree!" if(!$taxaInTree{$taxon});
    $tree->remove_Node($taxon);
  }
}

sub removeSingletons{
  my($tree,$remove,$settings)=@_;

  # A singleton is a node whose direct ancestor only has one child node.
  # Therefore for each leaf, find its direct ancestor and see if
  # it only has one child.  If so, I need to perform a splice.
  my @taxaInTree=$tree->get_leaf_nodes;
  for my $taxon(@taxaInTree){
    my $parentNode=$taxon->ancestor();
    my @siblingNode=$parentNode->get_all_Descendents();
    next if(@siblingNode > 1);

    # At this point, the node is an only child and so its parent
    # has to be spliced.
    $tree->splice(-remove=>[$parentNode],-preserve_lengths=>1);
  }
}

sub usage{
  local $0=basename($0);
  "$0: Removes a taxon from a tree
  Usage: $0 --tree tree.dnd taxon1 [taxon2...]
  "
}
