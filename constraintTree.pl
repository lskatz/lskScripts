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
  GetOptions($settings,qw(help tsv=s verbose reroot! min-taxa=i)) or die $!;
  $$settings{reroot}//=1;
  $$settings{'min-taxa'}||=1;
  my @tree=@ARGV;

  die usage() if($$settings{help});
  die "ERROR: need tsv" if(!$$settings{tsv});
  die "ERROR: need trees" if(!@tree);

  print join("\t",qw(file baseTaxon levelFromRoot numTaxa Sn Sp Score))."\n";
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
      my $m=constraintTree($treeObj,$$settings{tsv},$settings);
      if(!$$m{baseTaxon}){
        logmsg "WARNING: skipping tree in $t";
        next;
      }
      print join("\t",basename($t),$$m{baseTaxon},$$m{levelsFromRoot},scalar(@{$$m{sisterTaxa}}), $$m{Sn}, $$m{Sp}, $$m{Snsp})."\n";
    }
  }
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
    numTaxa=>0,
  );

  if($$settings{reroot}){
    eval{
      reroot($treeObj);
    };
    if($@){
      return \%return;
    }
  }
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

  # now that "outbreak" has been applied, calculate Sn or Sp
  my(%sn,%sp,%snsp);
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

      # Sensitivity and Specificity calculation
      $sn{$node->id}[$i]=$TP/($TP+$FN);
      $sp{$node->id}[$i]=$TN/($TN+$FP);
      $snsp{$node->id}[$i] = ($sn{$node->id}[$i]+$sp{$node->id}[$i])/2;

      # Add a penalty for each unknown in the outbreak clade
      $snsp{$node->id}[$i] = $snsp{$node->id}[$i] * ($knowns/($knowns+$unknowns));

      # Format a few values
      for($snsp{$node->id}[$i]){
        $_=sprintf("%0.2f",$_);
      }

      # I mean, why even bother if the score is super low.
      # Avoid a skewed score.
      # TODO avoid a skewed score with some kind of
      # exponential penalty.
      next if($sn{$node->id}[$i] < 0.01);
      next if($sp{$node->id}[$i] < 0.01);

      next if($numDescendents < $$settings{'min-taxa'});

      if($snsp{$node->id}[$i] > $return{Snsp} 
        || ($snsp{$node->id}[$i] == $return{Snsp} && $numDescendents > $return{numTaxa})
      ){
        $return{baseTaxon}=$node->id;
        $return{levelsFromRoot}=$i;
        $return{Snsp}=$snsp{$node->id}[$i];
        $return{mrca}=$ancestory[$i];
        $return{Sn}=$sn{$node->id}[$i];
        $return{Sp}=$sp{$node->id}[$i];
        $return{numTaxa}=$numDescendents;
        if($$settings{verbose}){
          logmsg "Taxon $return{baseTaxon} at level $return{levelsFromRoot} with score $return{Snsp} and $return{numTaxa} descendents. $unknowns unknowns in the clade.";
        }
      }
    }
  }
  
  # Find all taxon names in the "best" clade
  my @sisterTaxa=();
  if($return{mrca}){
    $return{sisterTaxa} = [map{$_->id} grep{$_->is_Leaf} $return{mrca}->get_Descendents];
  }

  #for(qw(Sn Sp Snsp)){
  #  $return{$_}//=0;
  #}
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
  The spreadsheet must be two columns: taxon-name and boolean (1 or 0)

  Usage: $0 --tsv file.tsv tree1.dnd [tree2.dnd...]
  --verbose
  --noreroot    Do not midpoint root
  --min-taxa  1 Number of taxa that must be in the target clade
  "
}
