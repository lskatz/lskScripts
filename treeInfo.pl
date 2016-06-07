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
  GetOptions($settings,qw(help taxa confidence|bootstrap all cutoff=i)) or die $!;
  $$settings{cutoff}||=0;

  # If the user wants all output then set it.
  if($$settings{all}){
    $$settings{$_}=1 for(qw(confidence taxa));
  }

  my @tree=@ARGV;

  die usage() if(!@tree || $$settings{help});

  my @header=qw(File numLeaves numNodes);
  push(@header,'confidence') if($$settings{confidence});
  push(@header,'taxa') if($$settings{taxa});
  print join("\t",@header)."\n";
  for my $t(@tree){
    treeInfo($t,$settings);
  }
  
  return 0;
}

sub treeInfo{
  my($tree,$settings)=@_;

  my $t=Bio::TreeIO->new(-file=>$tree)->next_tree;
  my @leaf=map{$_->id} $t->get_leaf_nodes;
  my @node=$t->get_nodes;


  my @out=($tree,scalar(@leaf),scalar(@node));
  if($$settings{confidence}){
    my $sum=0;
    my $count=0;
    my @bootstrap;
    for my $n(@node){
      next if($n->get_all_Descendents < 1);
      my $bootstrap=$n->bootstrap || $n->id;;
      next if(!defined($bootstrap));
      next if($bootstrap < $$settings{cutoff});
      push(@bootstrap,$bootstrap);
      $sum+=$bootstrap;
      $count++;
    }
    
    # does median or avg make more sense here?
    my $median=sprintf("%0.2f",median(@bootstrap));
    my $avg=sprintf("%0.2f",$sum/$count);
    push(@out,$avg);
  }

  push(@out,join(",",sort { $a cmp $b} @leaf)) if($$settings{taxa});

  print join("\t",@out)."\n";
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}


sub usage{
  $0=basename $0;
  "$0: gives information about a phylogeny
  Usage: $0 tree.dnd [tree2.dnd...]
  --taxa         Shows taxon names
  --confidence   Average confidence/bootstrap values
  --cutoff     0 Do not consider confidence values below this value
  --all          Shows all possible output fields
  "
}
