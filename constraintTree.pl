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
  GetOptions($settings,qw(help tree|trees|treeout=s tsv=s verbose reroot! min-taxa=i)) or die $!;
  $$settings{reroot}//=1;
  $$settings{'min-taxa'}||=1;
  $$settings{tree}||="";
  my @tree=@ARGV;

  die usage() if($$settings{help});
  die "ERROR: need tsv" if(!$$settings{tsv});
  die "ERROR: need trees" if(!@tree);

  my $treeout;
  if($$settings{tree}){
    $treeout=Bio::TreeIO->new(-file=>">$$settings{tree}");
  }

  print join("\t",qw(file baseTaxon levelsFromRoot numInClade Sn Sp Score outgroup))."\n";
  for my $t(@tree){
    if(!-e $t){
      die "ERROR: tree file doesn't exist: $t";
    }
    if(-s $t < 1){
      logmsg "Tree file is empty; skipping: $t";
      next;
    }
    eval{
      Bio::TreeIO->new(-file=>$t)->next_tree;
    };
    if($@){
      logmsg "Could not read file $t; skipping";
      next;
    }
    my $in=Bio::TreeIO->new(-file=>$t);
    while(my $treeObj=$in->next_tree){
      # Only look at the given root node...
      my @ancestorNodes = $treeObj->get_root_node;
      # ... unless the user supplies --reroot
      if($$settings{reroot}){
        @ancestorNodes=grep {!$_->is_Leaf} $treeObj->get_nodes;
      }

      my @allResults=();
      for my $node(@ancestorNodes){
        $treeObj->set_root_node($node);
        # Get the metrics of the tree into $m
        my $m=constraintTree($treeObj,$$settings{tsv},$settings);
        if(!$$m{baseTaxon}){
          next;
        }
        push(@allResults,$m);
      }
      if(!@allResults){
        logmsg "WARNING: skipping tree in $t";
        next;
      }

      # Sorted metrics array
      my @m = sort{$$b{Snsp} <=> $$a{Snsp} || $$b{Sn} <=> $$a{Sn} || $$b{Sp} <=> $$a{Sp}} @allResults;
      my %m=%{$m[0]}; # I don't feel like typing this whole variable name
      print join("\t",basename($t), $m{baseTaxon},$m{levelsFromRoot},scalar(@{$m{sisterTaxa}}), $m{Sn}, $m{Sp}, $m{Snsp});
      #print "\t".join(",",@{$m{outgroupTaxa}});
      print "\t".scalar(@{$m{outgroupTaxa}});
      print "\n";

      if($$settings{tree}){
        $treeout->write_tree($m{tree});
      }
    }
  }
  $treeout->close;

  # add an extra newline for readability
  open(my $fh, ">>", $$settings{tree}) or die $!;
  print $fh "\n";
  close $fh;
  
  return 0;
}


sub constraintTree{
  my($treeObj,$tsv,$settings)=@_;

  my %return=(
    baseTaxon=>"",
    levelsFromRoot=>0,
    Sn=>0,
    Sp=>0,
    Snsp=>0,
    sisterTaxa=>[],
    tree=>$treeObj,
    numTaxa=>0,
    outgroup=>treeOutgroup($treeObj,$settings),
  );
  $return{outgroupTaxa}=[sort{$a cmp $b} map{$_->id} grep{$_->is_Leaf} $return{outgroup}->get_all_Descendents];
  delete($return{outgroup}); # Not sure if I actually want to store this in the return hash

  my @node=();
  eval{
    @node=$treeObj->get_nodes;
  };
  if($@){
    return \%return;
  }

  my $inclusion=inclusionStatus($tsv,$settings);
  for my $node(@node){
    next if(!$node->is_Leaf);
    if(defined($$inclusion{$node->id})){
      #die "ERROR: can't find ".$node->id." in inclusion spreadsheet";
      $node->set_tag_value("outbreak",$$inclusion{$node->id});
    } else {
      $node->set_tag_value("outbreak",-1);
      logmsg "Warning: unknown inclusion status for ".$node->id;
    }
  }

  my $totalPositives=scalar(grep {$_->is_Leaf && $_->get_tag_values("outbreak")== 1} @node);
  my $totalNegatives=scalar(grep {$_->is_Leaf && $_->get_tag_values("outbreak")== 0} @node);
  my $totalUnknowns =scalar(grep {$_->is_Leaf && $_->get_tag_values("outbreak")==-1} @node);
  my $numNodes=@node;

  my @resultCombination; # array of hashes
  for my $node(@node){
    # For each leaf node, go up the ancestory chain to
    # record Sn and Sp
    next if(!$node->is_Leaf);

    my @ancestory=$treeObj->get_lineage_nodes($node);
    for(my $i=0;$i<@ancestory;$i++){
      # True positives, etc
      my @descendent=grep{$_->is_Leaf} $ancestory[$i]->get_Descendents;
      my $numDescendents=@descendent;
      my $TP = scalar(grep {$_->get_tag_values("outbreak")==1} @descendent);
      my $FP = scalar(grep {$_->get_tag_values("outbreak")==0} @descendent);
      my $TN = $totalNegatives - $FP;
      my $FN = $totalPositives - $TP;
      my $unknowns=scalar(grep {$_->get_tag_values("outbreak")==-1} @descendent);
      my $knowns  =scalar(grep {$_->get_tag_values("outbreak")!=-1} @descendent);

      # Can't simply have a clade of unknowns
      next if($knowns < 1);
      # If there is simply a zero value here, then there
      # is no redeeming this tree.
      next if($TN+$FP==0 || $TP+$FN==0);

      # Sensitivity and Specificity calculation
      my $sn    = $TP/($TP+$FN);
      my $sp    = $TN/($TN+$FP);
      my $snsp  = ($sn + $sp)/2;

      # Add a penalty for each unknown in the outbreak clade
      $snsp = $snsp * ($knowns/($knowns+$unknowns));

      # Format a few values
      for($sn,$sp,$snsp){
        $_=sprintf("%0.2f",$_);
      }

      # I mean, why even bother if the score is super low.
      # Avoid a skewed score.
      # TODO avoid a skewed score with some kind of
      # exponential penalty.
      next if($sn < 0.01);
      next if($sp < 0.01);

      next if($numDescendents < $$settings{'min-taxa'});

      my $treeCopy = Bio::Tree::Tree->new(-root=>$treeObj->get_root_node, -nodelete=>1);
      push(@resultCombination,{
        baseTaxon      => $node->id,
        levelsFromRoot => $i,
        Snsp           => $snsp,
        mrca           => $ancestory[$i],
        Sn             => $sn,
        Sp             => $sp,
        numTaxa        => $numDescendents,
        TP             => $TP,
        FP             => $FP,
        TN             => $TN,
        FN             => $FN,
        knowns         => $knowns,
        unknowns       => $unknowns,
        tree           => $treeCopy,
      });
    }
  }
  
  # If no results have been obtained, return the 
  # empty result.
  if(!@resultCombination){
    return \%return;
  }

  # Sort the results by the score, followed by sn, then
  # sp, and then finally the number of taxa in the 
  # target clade.
  my @sortedResults = sort{
    $$b{Snsp}    <=> $$a{Snsp} ||
    $$b{Sn}      <=> $$a{Sn}   ||
    $$b{Sp}      <=> $$a{Sp}   ||
    $$b{numTaxa} <=> $$a{numTaxa}
  } @resultCombination;

  for(qw(Snsp Sn Sp numTaxa baseTaxon levelsFromRoot mrca TP FP TN FN knowns unknowns tree)){
    $return{$_}=$sortedResults[0]{$_};
  }

  #logmsg "Taxon $return{baseTaxon} at level $return{levelsFromRoot} with score $return{Snsp} and $return{numTaxa} descendents. $return{unknowns} unknowns in the clade.";

  # Find all taxon names in the "best" clade
  my @sisterTaxa=();
  if($return{mrca}){
    $return{sisterTaxa} = [map{$_->id} grep{$_->is_Leaf} $return{mrca}->get_Descendents];
  }

  return \%return;
}

sub inclusionStatus{
  my($tsv,$settings)=@_;

  my %inclusion;
  open(my $fh, $tsv) or die "ERROR: could not read $tsv: $!";
  while(<$fh>){
    chomp;
    my ($isolateName,$bool)=split /\t/;
    # Force a boolean integer: the !! is "not-not" which
    # forces a 1 or empty value.  The plus 0 forces int.
    # Therefore, any "true" value evaluates to 1 and 
    # any "empty" value evaluates to 0.
    $inclusion{$isolateName} = !!$bool + 0;
  }
  close $fh;
  return \%inclusion;
}

#sub reroot{
#  my($tree,$settings)=@_;
#  
#  # Set a default outgroup before looking at the rest of the nodes.
#  my @node=$tree->get_nodes;
#  my $numNodes=@node;
#  my $outgroup=$node[0];
#  my $longest=$node[0]->branch_length || 0;
#
#  for(my $i=0;$i<$numNodes;$i++){
#    for(my $j=$i+1;$j<$numNodes;$j++){
#      if($tree->distance([$node[$i],$node[$j]]) > $longest){
#        $longest=$tree->distance([$node[$i],$node[$j]]);
#        $outgroup=$node[$i];
#      }
#    }
#  }
#  
#  $tree->reroot_at_midpoint($outgroup,'MidpointRoot');
#  return $outgroup;
#}

# Given the way the tree is rooted, return which node
# is the outgroup, defined by the fewest number of leaves
# of the first descendent of the root.
sub treeOutgroup{
  my($treeObj,$settings)=@_;

  # Get all direct descendents of the root node,
  # sorted by number of descendents.
  my @secondaryNode=sort{
    my @taxaA=sort {$a->id cmp $b->id} grep{$_->is_Leaf} $a->get_all_Descendents;
    my @taxaB=sort {$a->id cmp $b->id} grep{$_->is_Leaf} $b->get_all_Descendents;
    my $numA=scalar(@taxaA);
    my $numB=scalar(@taxaB);

    return 0 if($numA==0 && $numB==0);

    if($numA != $numB){
      return $numA <=> $numB;
    }

    # If they have the same number of taxa, go with
    # alphabetical order to keep a stable sort
    for(my $i=0;$i<$numA;$i++){
      if($taxaA[$i] ne $taxaB[$i]){
        return $taxaA[$i] cmp $taxaB[$i];
      }
    }
    
    die "INTERNAL ERROR: somehow this tree bifurcates into clades with equally named leaves???";
  } $treeObj->get_root_node->each_Descendent;

  for my $node(@secondaryNode){
    my $numLeaves=scalar(grep{$_->is_Leaf} $node->get_all_Descendents);
    if($numLeaves > 0){
      return $node;
    }
  }
  
  # TODO
  # If there are no outgroup leaves, then 
  # look at the other secondary nodes

  #die "ERROR: no secondary node has leaves";
  return $secondaryNode[0];
}

sub usage{
  "$0: find sensitivity and specificity of a tree or trees, 
  The spreadsheet must be two columns: taxon-name and boolean (1 or 0)

  Usage: $0 --tsv file.tsv tree1.dnd [tree2.dnd...]
  --verbose
  --noreroot            Do not try to root the tree every which way
  --min-taxa  1         Number of taxa that must be in the target clade
  --trees     file.dnd  Place final rerooted tree(s) in this file
  "
}

