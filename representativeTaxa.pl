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
    for(my $j=$i+1; $j<$numnodes; $j++){
      $distance{$taxon[$i]->id}{$taxon[$j]->id} = $tree->distance(-nodes=>[$taxon[$i],$taxon[$j]]);
      $distance{$taxon[$j]->id}{$taxon[$i]->id} = $distance{$taxon[$i]->id}{$taxon[$j]->id};
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

sub usage{
  "$0: Find representative taxa in each tree. Assumes one
  tree per tree file.

  Usage: $0 [options] tree.dnd [tree2.dnd...]

  --cluster-distance  0.001  The max distance between every
                             taxon in a cluster
  "
}
