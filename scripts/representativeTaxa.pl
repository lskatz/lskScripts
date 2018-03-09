#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use Bio::TreeIO;

local $0 = basename $0;
sub logmsg { print STDERR "$0: @_\n";}

exit main();

sub main{

  my $settings={};
  GetOptions($settings,qw(help cluster-distance|distance|max-distance=f)) or die $!;
  $$settings{'cluster-distance'}||=0.001;

  die usage() if(!@ARGV);

  for my $file(@ARGV){
    findRepresentatives($file,$settings);
  }
  
  return 0;
}

sub findRepresentatives{
  my($file,$settings)=@_;

  my $tree=Bio::TreeIO->new(-file=>$file)->next_tree;

  my @taxon = grep {$_->is_Leaf} $tree->get_nodes();
  my $numnodes=@taxon;

  # Find distances between all genomes
  logmsg "Finding distances between all taxa";
  my %distance;
  for(my $i=0;$i<$numnodes;$i++){
    print STDERR ".";
    my $taxonName1=$taxon[$i]->id;
    for(my $j=$i+1; $j<$numnodes; $j++){
      my $taxonName2=$taxon[$j]->id;
      my $distance=distanceBetweenTwoNodes($tree,$taxon[$i],$taxon[$j]);
      $distance{$taxonName1}{$taxonName2} = $distance;
      $distance{$taxonName2}{$taxonName1} = $distance;
    }
  }
  print STDERR "\n";


  # Cluster the taxa by distance
  # The index of %cluster is the representative genome,
  # and all other genomes have to be within X distance
  # of it.
  my %cluster;
  my $cluster_counter=0;
  for(my $i=0;$i<$numnodes;$i++){
    my $taxonName=$taxon[$i]->id;
    my $is_representative_taxon=1;
    for my $representative (keys(%cluster)) {
      if($distance{$taxonName}{$representative} < $$settings{'cluster-distance'}){
        push(@{ $cluster{$representative} }, $taxonName);
        $is_representative_taxon=0;
        last;
      }
    }

    if($is_representative_taxon){
      $cluster{$taxonName}=[$taxonName];
    }
  }

  for my $members(values(%cluster)){
    print join("\t",@$members)."\n";
  }

  logmsg "Found ".scalar(keys(%cluster))." clusters";

}

# http://cpansearch.perl.org/src/CJFIELDS/BioPerl-1.007002/Bio/Tree/TreeFunctionsI.pm
#   -> sub distance
# without error checking to speed it up
sub distanceBetweenTwoNodes{
  my($tree,$node1,$node2)=@_;

    my $lca = $tree->get_lca($node1,$node2);
    my $cumul_dist = 0;
    foreach my $current_node ($node1,$node2){
      do {
        $cumul_dist += $current_node->branch_length;

        $current_node = $current_node->ancestor || last;

      } while($current_node ne $lca);
    }

    return $cumul_dist;
}

sub usage{
  "$0: Find representative taxa in each tree. Assumes one
  tree per tree file.

  Usage: $0 [options] tree.dnd [tree2.dnd...]

  --cluster-distance  0.001  The max distance between every
                             taxon in a cluster
  "
}
