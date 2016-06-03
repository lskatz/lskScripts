#!/usr/bin/env perl
use warnings;
use strict;
use Bio::Tree::RandomFactory;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;
use Data::Dumper qw/Dumper/;
use List::Util qw/shuffle/;
use List::MoreUtils qw/uniq/;

local $0=basename $0;
sub logmsg { print STDERR "$0: @_\n"}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help numTrees=i)) or die $!;
  $$settings{numTrees}||=1;

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

  # Create as many random trees as requested
  my $factory=Bio::Tree::RandomFactory->new(-taxa=>\@sampleName,-maxcount=>$$settings{numTrees});
  while(my $tree=$factory->next_tree){
    # Alter the nodes to my liking
    for my $node($tree->get_nodes){
      # Get random branch lengths from the original trees
      my $newBranchLength=$branchLength[int(rand($numBranchLength))];
      $node->branch_length($newBranchLength);
      if($node->is_Leaf){
        
      } 
      # Give it a high bootstrap because we're not randomizing that
      else {
        $node->bootstrap(100);
        $node->id(100);
      }
    }
    print $tree->as_text("newick")."\n";
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
  "$0: uses real trees to create random trees
  Usage: $0 realtree.dnd [realtree2.dnd...]
  "
}
