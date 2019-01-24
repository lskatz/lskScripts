#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Bio::TreeIO;
use Bio::Tree::Statistics;
use Getopt::Long;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;

  if(!@ARGV){
    die usage();
  }

  print STDERR "Reading in files\n";
  my $guideTree; # will be the first tree
  my @bs_tree=();
  my $i;
  for my $file(@ARGV){
    if(!-e $file|| !-s $file){
      print STDERR "Not found or empty: $file ";
      next;
    }
    my $in=Bio::TreeIO->new(-file=>$file,-format=>"newick");
    #while(my $tree=next_tree_fast($in)){
    while(my $tree=$in->next_tree){
      if(!$guideTree){
        $guideTree = $tree;
        next;
      }
      push(@bs_tree, $tree);
      print STDERR ".";
    }
  }

  if(!$guideTree){
    die "ERROR: no guide tree found";
  }
  if(!@bs_tree){
    die "ERROR: no jack knife trees found";
  }

  print STDERR "\n";
  print STDERR "Combining jack knife files\n";
  my $biostat = Bio::Tree::Statistics->new();
  #my $bsTree=$biostat->assess_bootstrap(\@bs_tree,$guideTree);
  my $bsTree = assess_bootstrap($biostat, \@bs_tree, $guideTree);
  print STDERR "Reading internal nodes\n";
  for my $node($bsTree->get_nodes){
    print STDERR ".";
    next if($node->is_Leaf);

    if(!$node->id){
      $node->id($node->bootstrap);
    }
  }
  print STDERR "\n";
  print $bsTree->as_text("newick")."\n";

  return 0;
}

sub assess_bootstrap{
   my ($self,$bs_trees,$guide_tree) = @_;
   my @consensus;
 
   # internal nodes are defined by their children
 
   my (%lookup,%internal);
   my $i = 0;
   for my $tree ( $guide_tree, @$bs_trees ) {
       # Do this as a top down approach, can probably be
       # improved by caching internal node states, but not going
       # to worry about it right now.
 
       my @allnodes = $tree->get_nodes;
       my @internalnodes = grep { ! $_->is_Leaf } @allnodes;
       for my $node ( @internalnodes ) {
           my @tips = sort map { $_->id } 
                      grep { $_->is_Leaf() } $node->get_all_Descendents;
           my $id = "(".join(",", @tips).")";
           if( $i == 0 ) {
               $internal{$id} = $node->internal_id;
           } else { 
               $lookup{$id}++;
           }
       }
       $i++;
   }
   $i--; # do not count the guide tree in the denominator

   my @save;
   for my $l ( keys %lookup ) {
       if( defined $internal{$l} ) {#&& $lookup{$l} > $min_seen ) {
           my $intnode = $guide_tree->find_node(-internal_id => $internal{$l});
           $intnode->bootstrap(sprintf("%d",100 * $lookup{$l} / $i));
       }
   }
   return $guide_tree;
}

sub usage{
  "Usage: $0 guidetree.dnd jackknife.dnd [jackknife2.dnd...] > tree_with_confidence.nwk
  Where each jack knife tree can have multiple entries and the output tree
  will be a single entry with confidence values."
}
