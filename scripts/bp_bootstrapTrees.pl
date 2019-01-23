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
    my $in=Bio::TreeIO->new(-file=>$file,-format=>"newick");
    while(my $tree=next_tree_fast($in)){
    #while(my $tree=$in->next_tree){
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
    die "ERROR: no bootstrap trees found";
  }

  print STDERR "\n";
  print STDERR "Combining jack knife files\n";
  my $biostat = Bio::Tree::Statistics->new();
  my $bsTree=$biostat->assess_bootstrap(\@bs_tree,$guideTree);
  print STDERR "Reading internal nodes\n";
  for my $node($bsTree->get_nodes){
    print STDERR ".";
    next if($node->is_Leaf);

    #if(!$node->bootstrap){
    #  if($node->id){
    #    $node->bootstrap($node->id);
    #  }
    #}
    if(!$node->id){
      $node->id($node->bootstrap);
    }
  }
  print STDERR "\n";
  print $bsTree->as_text("newick")."\n";

  # explicitly end the script here to avoid confusion among all the subroutines
  exit 0;

  sub despace{
    my $dirty = shift; $dirty =~ s/\s+//gs; return $dirty;
  }
  sub dequote{
    my $dirty = shift;
    $dirty =~ s/^"?\s*(.+?)\s*"?$/$1/;
    return $dirty;
  }
  sub next_tree_fast {
      my ($treeio) = @_;
      local $/ = ";\n";
      return unless $_ = $treeio->_readline;

      s/[\r\n]//gs;

      my $score;
      s/([^"]*)(".+?")([^"]*)/despace->($1) . dequote->($2) . despace->($3)/egsx;
   
      if (s/^\s*\[([^\]]+)\]//) {
          my $match = $1;
          $match =~ s/\s//g;
          $match =~ s/lh\=//;
          if ( $match =~ /([-\d\.+]+)/ ) {
              $score = $1;
          }
      }
   
      $treeio->_eventHandler->start_document;
      # Call the parse_newick method as defined in NewickParser.pm
      $treeio->parse_newick($_);
      my $tree = $treeio->_eventHandler->end_document;
     
      # Add the tree score afterwards if it exists.
      if (defined $tree) {
        $tree->score($score);
        return $tree;
      }
      return $tree;
  }

  return 0;
}

sub usage{
  "Usage: $0 guidetree.dnd jackknife.dnd [jackknife2.dnd...] > bootstraps.dnd
  Where each jack knife tree can have multiple entries and the output tree
  will be a single entry with bootstraps."
}
