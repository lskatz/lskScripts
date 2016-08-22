#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;
use List::MoreUtils qw/uniq/;

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help taxa confidence|bootstrap all cutoff=i diff)) or die $!;
  $$settings{cutoff}||=0;

  # If the user wants all output then set it.
  if($$settings{all}){
    $$settings{$_}=1 for(qw(confidence taxa));
  }

  my @tree=@ARGV;
  my $numTrees=@tree;

  die usage() if(!@tree || $$settings{help});

  my @header=qw(File numLeaves numNodes);
  push(@header,'confidence') if($$settings{confidence});
  push(@header,'taxa') if($$settings{taxa});
  print join("\t",@header)."\n";

  my %taxa;
  my @taxa;
  for my $t(@tree){
    my $taxa=treeInfo($t,$settings);
    for(@$taxa){
      push(@taxa,@$taxa);
      $taxa{$t}{$_}++;
      if($taxa{$t}{$_} > 1){
        logmsg "WARNING: taxon $_ was found more than once in $t";
      }
    }
  }
  @taxa=uniq(sort {$a cmp $b } @taxa);

  # Do any of the trees show differences?
  if($$settings{diff}){
    my %refTaxa=%{ $taxa{$tree[0]} };
    my $refNumTaxa=scalar(keys(%refTaxa));

    #TODO figure out a better way for multiple comparisons
    # but for now we will just do two trees.
    #my $numTrees=2;
    my %diff;

    # header
    print "\nDIFF\n";
    print join("\t",@tree)."\n";
    my $numTreeHeader="";
    $numTreeHeader.="numTaxa=".scalar( keys(%{$taxa{$_}}) ) ."\t" for(@tree);
    $numTreeHeader=~s/\t$//;
    print $numTreeHeader."\n";

    for(my $i=1;$i<$numTrees;$i++){
      my @queryTaxa=keys(%{$taxa{$tree[$i]}});

      # Find extra taxa
      for my $taxon(@queryTaxa){
        if(!$refTaxa{$taxon}){
          $diff{$taxon}{$tree[$i]}=$taxon;
          $diff{$taxon}{$tree[0]}="*";
        }
      }

      # Find missing taxa
      my %queryTaxa;
      @queryTaxa{@queryTaxa}=(1) x scalar(@queryTaxa);
      for my $refTaxon(keys(%refTaxa)){
        if(!$queryTaxa{$refTaxon}){
          $diff{$refTaxon}{$tree[$i]}="*";
          $diff{$refTaxon}{$tree[0]}=$refTaxon;
        }
      }
    }

    for my $taxon(@taxa){
      my $diffLine="";
      for my $tree(@tree){
        $diff{$taxon}{$tree}//=$taxon;
        $diffLine.=$diff{$taxon}{$tree}."\t";
      }
      $diffLine=~s/\t$//;
      print $diffLine."\n";
    }

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
  return \@leaf;
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
  --diff         Shows information about differences between the
                 first tree and the others
  "
}
