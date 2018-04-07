#!/usr/local/bin/env perl
use strict;
use warnings;
use Bio::Perl;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;
use Math::Round qw/ceil/;

local $0=basename($0);
sub logmsg{print STDERR "$0: @_\n"}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help));

  my $$settings{bootstrap}//=70;
  die usage() if($$settings{help});

  my @tree=@ARGV;
  for my $tree(@tree){
    my $in=Bio::TreeIO->new(-file=>$tree);
    while(my $treeObj=$in->next_tree){
      $treeObj=splitPolytomies($treeObj,undef,$settings);
    }
  }
  
  return 0;
}

sub splitPolytomies{
  my($treeObj,$node,$settings)=@_;

  my $node ||= $treeObj->get_root_node;

  my @descs = $node->each_Descendent;
  if (@descs > 2) {
    # Many nodes have no identifying names, a simple warning is probably
    # enough.

    $treeObj->warn("Node has more than two descendants\nWill do an arbitrary balanced split");
    my @working = @descs;
    # create an even set of artifical nodes on which to later hang the descs
    my $half = ceil(@working / 2);
    my @artificials;
    while ($half > 1) {
      my @this_level;
      foreach my $top_node (@artificials || $node) {
        for (1..2) {
          my $art = $top_node->new(-id => "artificial_".++$treeObj->{_art_num});
          $top_node->add_Descendent($art);
          push(@this_level, $art);
        }
      }
      @artificials = @this_level;
      $half--;
    }
    # attach two descs to each artifical leaf
    foreach my $art (@artificials) {
      for (1..2) {
        my $desc = shift(@working) || $node->new(-id => "artificial_".++$treeObj->{_art_num});
        $desc->ancestor($art);
      }
    }
  }
  elsif (@descs == 1) {
    # ensure that all nodes have 2 descs
    $node->add_Descendent($node->new(-id => "artificial_".++$treeObj->{_art_num}));
  }
  # recurse
  foreach my $desc (@descs) {
    splitPolytomies($treeObj,$desc,$settings);
  }

  return $treeObj;
}

sub usage{
  "$0: split polytomies in a predictable way
  Usage: $0 tree.dnd [tree2.dnd...] > out.dnd
  --bootstrap  70  The minimum bootstrap value where a
                   clade will be considered a polytomy
  "
}
