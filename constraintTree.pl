#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw/basename/;
use Getopt::Long qw/GetOptions/;
use Data::Dumper;

use Bio::TreeIO;

local $0=basename $0;

sub logmsg{print STDERR "$0: @_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tsv=s)) or die $!;
  my @tree=@ARGV;

  die usage() if($$settings{help});
  die "ERROR: need tsv" if(!$$settings{tsv});
  die "ERROR: need trees" if(!@tree);

  for my $t(@tree){
    my $in=Bio::TreeIO->new(-file=>$t);
    while(my $treeObj=$in->next_tree){
      my $metrics=constraintTree($treeObj,$$settings{tsv},$settings);
      die Dumper $metrics;
    }
  }
}

sub constraintTree{
  my($treeObj,$tsv,$settings)=@_;

  reroot($treeObj);

  my $inclusion=inclusionStatus($tsv,$settings);
  for my $node($treeObj->get_nodes){
    next if(!$node->is_Leaf);
    $node->set_tag_value("outbreak",$$inclusion{$node->id});
  }

  my $totalPositives=scalar(grep {$_->is_Leaf && $_->get_tag_values("outbreak")==1} $treeObj->get_nodes);
  my $totalNegatives=scalar(grep {$_->is_Leaf && $_->get_tag_values("outbreak")==0} $treeObj->get_nodes);

  # now that "outbreak" has been applied, see Sn or Sp
  my(%sn,%sp,%snsp);
  # Also try to find the best Sn/Sp along the way
  my $winningTaxon="";
  my $winningLevel=0;
  my $winningSnsp=0;
  for my $node($treeObj->get_nodes){
    # For each leaf node, go up the ancestory chain to
    # record Sn and Sp
    next if(!$node->is_Leaf);

    my @ancestory=$treeObj->get_lineage_nodes($node);
    for(my $i=0;$i<@ancestory;$i++){
      # True positives, etc
      my $TP = scalar(grep {$_->is_Leaf && $_->get_tag_values("outbreak")==1} $ancestory[$i]->get_Descendents);
      my $FP = scalar(grep {$_->is_Leaf && $_->get_tag_values("outbreak")==0} $ancestory[$i]->get_Descendents);
      my $TN = $totalNegatives - $FP;
      my $FN = $totalPositives - $TP;

      # Sensitivity and Specificity calculation
      $sn{$node->id}[$i]=$TP/($TP+$FN);
      $sp{$node->id}[$i]=$TN/($TN+$FP);
      $snsp{$node->id}[$i] = ($sn{$node->id}[$i]+$sp{$node->id}[$i])/2;

      # I mean, why even bother if the score is super low
      next if($sn{$node->id}[$i] < 0.01);
      next if($sp{$node->id}[$i] < 0.01);

      if($snsp{$node->id}[$i] > $winningSnsp){
        $winningTaxon=$node->id;
        $winningLevel=$i;
        $winningSnsp=$snsp{$node->id}[$i];
        logmsg $winningTaxon,$winningLevel,$winningSnsp,$TP,$FP,$TN,$FN;
      }
    }
  }

  return {baseTaxon=>$winningTaxon, levelsFromRoot=>$winningLevel, Sn=>$sn{$winningTaxon}[$winningLevel], Sp=>$sp{$winningTaxon}[$winningLevel], Snsp=>$snsp{$winningTaxon}[$winningLevel]};
}

sub inclusionStatus{
  my($tsv,$settings)=@_;

  my %inclusion;
  open(my $fh, $tsv) or die "ERROR: could not read $tsv: $!";
  while(<$fh>){
    chomp;
    my ($dataset,$event,$isolateName,$sra_acc)=split /\t/;
    $inclusion{$isolateName} = !!$event + 0; # force int
    #if($event){
    #  $inclusion{$isolateName}="true";
    #} else {
    #  $inclusion{$isolateName}="false";
    #}
  }
  close $fh;
  return \%inclusion;
}

sub reroot{
  my($tree,$settings)=@_;
  
  # Set a default outgroup before looking at the rest of the nodes.
  my @node=$tree->get_nodes;
  my $numNodes=@node;
  my $outgroup=$node[0];
  my $longest=$node[0]->branch_length || 0;

  for(my $i=0;$i<$numNodes;$i++){
    for(my $j=$i+1;$j<$numNodes;$j++){
      if($tree->distance([$node[$i],$node[$j]]) > $longest){
        $longest=$tree->distance([$node[$i],$node[$j]]);
        $outgroup=$node[$i];
      }
    }
  }
  $tree->reroot_at_midpoint($outgroup,'MidpointRoot');
  return $outgroup;
}

sub usage{
  "$0: find sensitivity and specificity of a tree or trees, 
  given a spreadsheet of inclusion/exclusion information

  Usage: $0 --tsv file.tsv tree1.dnd [tree2.dnd...]
  "
}
