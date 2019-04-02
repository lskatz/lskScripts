#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use Getopt::Long qw/GetOptions/;
use Data::Dumper qw/Dumper/;
use Scalar::Util qw/looks_like_number/;

use Bio::Tree::Draw::Cladogram;
use Bio::TreeIO;

local $0 = basename $0;
sub logmsg {print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings, qw(help tempdir outfile=s)) or die $!;
  $$settings{tempdir} ||= tempdir("$0.XXXXXX", TMPDIR=>1, CLEANUP=>1);

  my($tree1, $tree2) = @ARGV;

  my $t1 = readTree($tree1,$settings);
  my $t2 = readTree($tree2,$settings);

  logmsg "Creating tanglegram";
  my $tanglegram = Bio::Tree::Draw::Cladogram->new(-tree=>$t1, -second=>$t2, -compact=>0, -ratio=>1, -column=>200);
  
  # Print to temp file
  my $eps = "$$settings{tempdir}/tanglegram.eps";
  logmsg "Creating eps representation at $eps";
  $tanglegram->print(-file=>$eps);
  # Correct the bioperl module output
  open(my $fh1, $eps) or die "ERROR: could not read $eps: $!";
  open(my $fh2, ">", "$eps.eps") or die "ERROR: could not write to $eps.eps: $!";
  my($prevX,$prevY)=(0,0);
  while(<$fh1>){
    chomp;
    my @F=split /\s+/;
    if(/(lineto|moveto)/){
      if(looks_like_number($F[0])){
        $prevX=$F[0];
      } else {
        $F[0]=$prevX;
      } if(looks_like_number($F[1])){
        $prevY=$F[1];
      } else {
        $F[1]=$prevY;
      }
    } 
    print $fh2 join(" ", @F)."\n";
  }
  close $fh1;
  close $fh2;

  $eps = "$eps.eps";

  if($$settings{outfile}){
    logmsg "Printing to $$settings{outfile}";
    system("convert $eps $$settings{outfile}");
    die "ERROR: could not convert $eps to $$settings{outfile}: $!" if $?;
  } else {
    open(my $fh, $eps) or die "ERROR: could not read $eps: $!";
    while(<$fh>){
      print;
    }
    close $fh;
  }

  return 0;
}

sub readTree{
  my($t,$settings)=@_;
  my $obj = Bio::TreeIO->new(-file=>$t)->next_tree;
  
  # Reroot at longest branch length
  my $longestBranchLength = 0;
  my $newRootNode;
  for my $node($obj->get_nodes){
    if(defined($node->branch_length) && $node->branch_length > $longestBranchLength){
      $longestBranchLength = $node->branch_length;
      $newRootNode = $node;
    }
  }
  $obj->reroot_at_midpoint($newRootNode);

  return $obj;
}

sub usage{
  "Draws a tanglegram between two trees
  Usage: $0 [options] t1.dnd t2.dnd > out.eps
  --outfile  ''  If specified, the image will not go to
                 stdout.  The extension will dictate
                 the file format.
  "
}
