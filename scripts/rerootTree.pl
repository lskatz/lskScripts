#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;
use List::Util qw/sum/;
use Statistics::Descriptive;
use File::Temp qw/tempdir tempfile/;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help root-on|root-with|root=s)) or die $!;

  $$settings{'root-on'} || die "ERROR: parameter --root-on is required. --help for more information.";

  die usage() if(!@ARGV || $$settings{help});

  my(@query)=@ARGV;
  my $treeout=Bio::TreeIO->new(-format=>"newick");

  # For each query, reroot and print
  for my $q(@query){
    my $treein = Bio::TreeIO->new(-file=>$q,-format=>"newick");
    my $treeCounter=0;
    while(my $treeObject = $treein->next_tree){
      $treeCounter++;
      my @leaves = sort {$a->id cmp $b->id} grep {$_->is_Leaf} $treeObject->get_nodes;
      if(scalar(@leaves) < 2){
        if(scalar($treeObject->get_nodes) < 3){
          logmsg "Skipping: only ".scalar($treeObject->get_nodes)." nodes found in tree $treeCounter in $q";
          logmsg "  Possible reason: empty line in tree file";
          next;
        }
        die "ERROR: there are fewer than 2 leaves on tree $treeCounter in $q:\n".Dumper [map{$_->id} @leaves];
      }
      my @node = grep {$_->id eq $$settings{'root-on'}} @leaves;
      if(@node > 1){
        die "ERROR: found multiple nodes named ".$$settings{'root-on'}." in $q";
      }
      logmsg "Rerooting tree $treeCounter in $q";
      
      my $was_rerooted=$treeObject->reroot($node[0]);
      if(!$was_rerooted){
        die "ERROR: could not reroot tree $treeCounter in $q";
      }

      $treeout->write_tree($treeObject);

    }
  }

  return 0;
}

sub usage{
  "$0: Roots a set of trees on the same leaf.
  Output trees will be in the same order as tree parameters
  Usage: $0 --root-on LEAF tree.dnd [tree2.dnd...] > trees.dnd
  --root-on    ''  The name of the leaf
  "
}

