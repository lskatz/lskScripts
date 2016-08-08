#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use File::Copy qw/cp mv/;
use File::Spec;
use Getopt::Long;
use Statistics::Descriptive;
use Math::Gauss qw/cdf pdf/;

sub logmsg{local $0=basename $0; print STDERR "$0 @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s method=s numcpus=i numtrees=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{numtrees}||=1000;
  $$settings{method}||="kf";
  $$settings{tempdir}||=tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{tempdir}=File::Spec->rel2abs($$settings{tempdir});

  my($ref,@query)=@ARGV;

  die usage() if(!$ref || !@query || $$settings{help});

  # Check executables
  for my $exe (qw(treedist randTrees.pl)){
    system("which $exe >/dev/null 2>&1");
    die "ERROR: I did not find $exe in your PATH" if $?;
  }

  logmsg "Temporary dir is $$settings{tempdir}";

  print join("\t",qw(Ref Query num obs avg stdev Z p))."\n";
  for my $query(@query){
    my $stat={};
    $$settings{method}=lc($$settings{method});
    if($$settings{method} eq "kf"){
      $stat=kuhnerFelsenstein($ref,$query,$settings);
    } elsif($$settings{method} eq "rf"){
      $stat=robinsonFoulds($ref,$query,$settings);
    } else {
      die "ERROR: I do not understand method $$settings{method}";
    }

    $$stat{$_}=sprintf("%0.2f",$$stat{$_}) for(qw(obs avg stdev Z p));
    print join("\t",$$stat{ref},$$stat{query}, $$stat{num}, $$stat{obs},
                    $$stat{avg}, $$stat{stdev}, $$stat{Z}, $$stat{p});
    print "\n";
  }

  return 0;
}

sub kuhnerFelsenstein{
  my($ref,$query,$settings)=@_;

  # Clean out certain files in the tempdir that may
  # or may not be there.
  # Unfortunately treedist has problems with overwriting
  # files.
  for("$$settings{tempdir}/intree","$$settings{tempdir}/intree2","$$settings{tempdir}/outfile"){
    unlink($_);
  }
  
  cp($ref,$$settings{tempdir}) or die $!;
  cp($query,$$settings{tempdir}) or die $!;

  my $tmpRefTree="$$settings{tempdir}/".basename($ref);
  my $tmpQueryTree="$$settings{tempdir}/".basename($query);

  
  # Create a list of trees to compare against in $comparisonTrees:
  #   The first tree is the reference tree
  #   The next trees are randomly made trees from 
  #   the query tree.
  #   Therefore, the question being answered is, is the query
  #   closer to the ref than random trees?
  my $comparisonTrees="$$settings{tempdir}/compareAgainst.dnd";
  cp($tmpRefTree,$comparisonTrees) or die $!;
  logmsg "Generating $$settings{numtrees} random trees";
  system("randTrees.pl $tmpQueryTree --force-binary --numcpus $$settings{numcpus} --numtrees $$settings{numtrees} >> $comparisonTrees"); die if $?;
  # We need the files to be named intree and intree2 because of inflexible treedist
  symlink(basename($tmpQueryTree),"$$settings{tempdir}/intree");
  symlink(basename($comparisonTrees),"$$settings{tempdir}/intree2");


  system("cd $$settings{tempdir} && echo -ne '2\nL\nS\nY\n' | treedist > treedist.log 2> treedist.out");
  die "ERROR with treedist ".`cat $$settings{tempdir}/treedist.log` if $?;

  # Read the output file.
  # "1" in the outfile indicates the query tree.
  # "2" in the outfile indicates the ref tree.
  # >2 indicates a random tree.
  # I only care about query vs ref and query vs randoms.
  # I do not care about anything vs query.
  my @score;  # background distribution
  my $obs;    # observed value
  open(TREEDISTOUT,"$$settings{tempdir}/outfile") or die "ERROR: could not read treedist output: $!";
  while(<TREEDISTOUT>){
    chomp;
    my($tree1,$tree2,$score)=split(/\s+/,$_);
    next if($tree1 ne 1);
    
    if($tree2 eq 2){
      $obs=$score;
      next;
    }
    push(@score,$score);
  }
  close TREEDISTOUT;
  
  my $stat=Statistics::Descriptive::Full->new();
  $stat->add_data(@score);
  my $avg=$stat->mean;
  my $stdev=$stat->standard_deviation;
  
  my $Z=($obs - $avg)/$stdev;
  my $p=cdf($obs,$avg,$stdev);

  return {Z=>$Z, p=>$p, obs=>$obs, num=>scalar(@score), avg=>$avg, stdev=>$stdev, query=>$query, ref=>$ref};
}

sub robinsonFoulds{
  my($ref,$query,$settings)=@_;

  # Clean out certain files in the tempdir that may
  # or may not be there.
  # Unfortunately treedist has problems with overwriting
  # files.
  for("$$settings{tempdir}/intree","$$settings{tempdir}/intree2","$$settings{tempdir}/outfile"){
    unlink($_);
  }
  
  cp($ref,$$settings{tempdir}) or die $!;
  cp($query,$$settings{tempdir}) or die $!;

  my $tmpRefTree="$$settings{tempdir}/".basename($ref);
  my $tmpQueryTree="$$settings{tempdir}/".basename($query);

  
  # Create a list of trees to compare against in $comparisonTrees:
  #   The first tree is the reference tree
  #   The next trees are randomly made trees from 
  #   the query tree.
  #   Therefore, the question being answered is, is the query
  #   closer to the ref than random trees?
  my $comparisonTrees="$$settings{tempdir}/compareAgainst.dnd";
  cp($tmpRefTree,$comparisonTrees) or die $!;
  logmsg "Generating $$settings{numtrees} random trees";
  system("randTrees.pl $tmpQueryTree --numcpus $$settings{numcpus} --numtrees $$settings{numtrees} >> $comparisonTrees"); die if $?;
  # We need the files to be named intree and intree2 because of inflexible treedist
  symlink(basename($tmpQueryTree),"$$settings{tempdir}/intree");
  symlink(basename($comparisonTrees),"$$settings{tempdir}/intree2");


  system("cd $$settings{tempdir} && echo -ne 'D\n2\nL\nS\nY\n' | treedist > treedist.log 2> treedist.out");
  die "ERROR with treedist ".`cat $$settings{tempdir}/treedist.log` if $?;

  # Read the output file.
  # "1" in the outfile indicates the query tree.
  # "2" in the outfile indicates the ref tree.
  # >2 indicates a random tree.
  # I only care about query vs ref and query vs randoms.
  # I do not care about anything vs query.
  my @score;  # background distribution
  my $obs;    # observed value
  open(TREEDISTOUT,"$$settings{tempdir}/outfile") or die "ERROR: could not read treedist output: $!";
  while(<TREEDISTOUT>){
    chomp;
    my($tree1,$tree2,$score)=split(/\s+/,$_);
    next if($tree1 ne 1);
    
    if($tree2 eq 2){
      $obs=$score;
      next;
    }
    push(@score,$score);
  }
  close TREEDISTOUT;
  
  my $stat=Statistics::Descriptive::Full->new();
  $stat->add_data(@score);
  my $avg=$stat->mean;
  my $stdev=$stat->standard_deviation;
  
  my $Z=($obs - $avg)/$stdev;
  my $p=cdf($obs,$avg,$stdev);

  return {Z=>$Z, p=>$p, obs=>$obs, num=>scalar(@score), avg=>$avg, stdev=>$stdev, query=>$query, ref=>$ref};
}

sub usage{
  local $0=basename $0;
  "$0: compares two trees
  Usage: $0 ref.dnd query.dnd
  --method    kf     kf: kuhner-felsenstein
                     rf: robinson-foulds
  --numtrees  1000   How many random trees to compare against
  "
}

