#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename qw/basename/;
use List::Util qw/sum/;
use List::MoreUtils qw/uniq/;

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help list=s taxa confidence|bootstrap min-confidence all cutoff=f diff)) or die $!;
  $$settings{cutoff}||=0;

  # If the user wants all output then set it.
  if($$settings{all}){
    $$settings{$_}=1 for(qw(confidence taxa));
  }

  my @tree=@ARGV;
  my $numTrees=@tree;

  die usage() if(!@tree || $$settings{help});

  my @header=qw(File numLeaves numNodes avgBranchLength);
  push(@header,'confidence') if($$settings{confidence});
  push(@header,'min-confidence') if($$settings{'min-confidence'});
  push(@header,'taxa') if($$settings{taxa});
  print join("\t",@header)."\n";

  my %taxa;
  my @taxa;
  for my $t(@tree){
    if(my $thing_to_list = $$settings{list}){
      my $listing = treeCategoryList($thing_to_list, $t, $settings);
      for my $item(@$listing){
        print $item."\n";
      }
      next;
    }
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
  my $unnamedNodeCounter=0;
  my @leaf=map{$_->id || "unnamedTaxon".++$unnamedNodeCounter} $t->get_leaf_nodes;
  my @node=$t->get_nodes;

  my $avgBranchLength = sum(map{$_->branch_length||0} @node)/scalar(@node);

  my @out=($tree,scalar(@leaf),scalar(@node),$avgBranchLength);
  if($$settings{confidence}){
    my $minConfidence=100;
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

      if($minConfidence > $bootstrap){
        $minConfidence=$bootstrap;
      }
    }
    
    # does median or avg make more sense here?
    my($median,$avg)=(0,0);
    if($count > 1){
      $median=sprintf("%0.2f",median(@bootstrap));
      $avg=sprintf("%0.2f",$sum/$count);
    }
    push(@out,$avg);

    if($$settings{'min-confidence'}){
      push(@out,$minConfidence);
    }
  }


  push(@out,join(",",sort { $a cmp $b} @leaf)) if($$settings{taxa});

  print join("\t",@out)."\n";
  return \@leaf;
}

sub treeCategoryList{
  my($thing_to_list, $tree, $settings) = @_;

  $thing_to_list = lc($thing_to_list);
  my @list;
  my $sortOrder = 'abc';

  my $t=Bio::TreeIO->new(-file=>$tree)->next_tree;
  # taxon names, bootstrap values, branch lengths
  for my $node($t->get_nodes){
    if($thing_to_list =~ /bootstrap|confidence/){
      next if($node->is_Leaf);
      push(@list, $node->bootstrap || $node->id || 100);
      $sortOrder = '123';
    } elsif($thing_to_list =~ /tax/){
      next if(! $node->is_Leaf);
      my $id = $node->id || '';
      push(@list, $id);
      $sortOrder = 'abc';
    } elsif($thing_to_list =~ /branch|length/){
      my $length = $node->branch_length || 0;
      push(@list, $length);
      $sortOrder = '123';
    } else {
      die "ERROR: you provided --list but I do not understand the value $thing_to_list";
    }
  }

  if($sortOrder eq '123'){
    @list = sort{$a <=> $b} @list;
  } else {
    @list = sort{$a cmp $b} @list;
  }

  return \@list;
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    elsif($len==0){
        return 0;
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
  --list      '' List either: taxa, bootstrap, branches
  --taxa         Shows taxon names
  --confidence   Average confidence/bootstrap values
  --cutoff     0 Do not consider confidence values below this value
  --min-confidence  Display the minimum confidence value
  --all          Shows all possible output fields
  --diff         Shows information about differences between the
                 first tree and the others
  "
}
