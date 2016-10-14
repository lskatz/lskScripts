#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use File::Temp qw/tempdir tempfile/;
use File::Copy qw/cp mv/;
use File::Spec;
use Getopt::Long;
use Statistics::Descriptive;
use Math::Gauss qw/cdf pdf/;

use Bio::Phylo::IO;
use Bio::Phylo::Factory;
use Bio::TreeIO;

use threads;

my $startTime=time();
sub logmsg{
  my $duration=time() - $startTime;
  local $0=basename $0; 
  my $TID='TID'.threads->tid;
  #print STDERR "$0 [$duration]: @_\n";
  print STDERR "$0($TID): @_\n";
}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s method=s numcpus=i numtrees=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{numtrees}||=0;
  $$settings{method}||="kf";
  $$settings{method}||=lc($$settings{method});
  $$settings{tempdir}||=tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{tempdir}=File::Spec->rel2abs($$settings{tempdir});

  my($ref,@query)=@ARGV;

  die usage() if(!$ref || !@query || $$settings{help});
  for($ref,@query){
    die "ERROR: cannot find $_" if(!-e $_);
  }

  # Check executables
  for my $exe (qw(randTrees.pl)){
    system("which $exe >/dev/null 2>&1");
    die "ERROR: I did not find $exe in your PATH" if $?;
  }

  logmsg "Temporary dir is $$settings{tempdir}";

  print join("\t",qw(Ref Query num obs avg stdev Z p))."\n";
  for my $query(@query){
    my $stat=bioPhyloDist($query,$ref,$$settings{method},$settings);

    # Round to two decimal places for most things
    $$stat{$_}=sprintf("%0.4f",$$stat{$_}) for(qw(obs avg stdev Z));
    # Scientific notation for p values
    $$stat{$_}=sprintf("%0.2e",$$stat{$_}) for(qw(p));
    print join("\t",$$stat{ref},$$stat{query}, $$stat{num}, $$stat{obs},
                    $$stat{avg}, $$stat{stdev}, $$stat{Z}, $$stat{p});
    print "\n";
  }

  return 0;
}

sub bioPhyloDist{
  my($query,$ref,$method,$settings)=@_;

  # First find the observed value.
  # The tree might be multifurcating and therefore might give a non-zero
  # observed result.  Rerun this test until we see a zero.
  my $obs;
  my $tries=0;
  do{
    $obs=observedScore($query,$ref,$method,$settings);
    if($tries++ > 10000){
      die "ERROR: tried to get a zero self vs self score for $ref but couldn't after $tries tries. The tree might be too multifurcating. Force bifurcation and try again";
    }
  } while($ref eq $query && $obs != 0);

  # Next run a test against a random background distribution

  my (@scores); # scores from the statistical test
  my ($avg,$stdev,$Z,$p)=(0) x 4;

  # Need more than 1 for a statistical test against
  # a distribution
  if($$settings{numtrees} > 1){
    # multithreaded random trees
    my @thr;
    my $treesPerThread=int($$settings{numtrees}/$$settings{numcpus})+1;
    for(0..$$settings{numcpus}-1){
      $thr[$_]=threads->new(\&randScores,$query,$ref,$$settings{method},$treesPerThread,$settings);
    }

    # Combine all random query trees. Results are found
    # in $randTreesQ
    for(@thr){
      my $scoresList=$_->join;
      # There was a problem with $scoresList and I put in this
      # debugging statement, but there really isn't a 
      # reason to remove it and so I kept it.
      die "INTERNAL ERROR: scoresList is not a list!" if(ref($scoresList) ne "ARRAY");
      push(@scores,@$scoresList);
    }
    # Respect the number of trees requested
    splice(@scores,$$settings{numtrees});

    my $stat=Statistics::Descriptive::Full->new();
    $stat->add_data(@scores);
    $avg=$stat->mean;
    $stdev=$stat->standard_deviation;
    $stdev=1e-9 if($stdev <= 0);
    $Z=($obs - $avg)/$stdev;
    $p=cdf($obs,$avg,$stdev);
  }
  
  return {Z=>$Z, p=>$p, obs=>$obs, num=>scalar(@scores), avg=>$avg, stdev=>$stdev, query=>$query, ref=>$ref};

}

sub randScores{
  my($query,$ref,$method,$numTrees,$settings)=@_;

  my $refObj=Bio::Phylo::IO->parse(
    -format=>"newick",
    -file=>$ref,
  )->first;

  logmsg "Making random trees from $query";
  my @newickString=`randTrees.pl $query --force-binary --numcpus 1 --numtrees $numTrees 2>/dev/null`;
  die "ERROR with randTrees.pl: $!" if $?;

  logmsg "Comparing random trees against $ref";
  my @dist;
  for my $newick(@newickString){
    my ($fh, $randQueryTree) = tempfile("randQuery.XXXXXX",DIR=>$$settings{tempdir},SUFFIX=>".dnd");
    print $fh $newick;
    close $fh;
    push(@dist,observedScore($randQueryTree,$ref,$method,$settings));
  }

  return \@dist;
}

sub observedScore{
  my($query,$ref,$method,$settings)=@_;

  my $queryObj=Bio::Phylo::IO->parse(
    -format=>"newick",
    -file=>$query,
  )->first;
  my $refObj=Bio::Phylo::IO->parse(
    -format=>"newick",
    -file=>$ref,
  )->first;

  for($queryObj,$refObj){
    # Break polytomies and make the tree binary.
    $_=$_->resolve()->deroot();

    # midpoint root: makes the score not match Phylip
    #my $midpointNode=$_->get_midpoint();
    #$_->reroot($midpointNode);
    #$_=$_->deroot();
  }

  # Set the tree objects equal to each other if the 
  # file paths are the same. Therefore the resolve()
  # function would result in the same tree.
  if($query eq $ref){
    $queryObj=$refObj;
  }

  #logmsg "$refObj $queryObj";
  #if(!$refObj->is_binary()){
  #  die "ERROR: the ref is not binary!\n  file: $ref\n  ".`cat $ref`;
  #}
  #if(!$queryObj->is_binary()){
  #  die "ERROR: the query is not binary!\n  file: $query\n  ".`cat $query`;
  #}

  # Find the score based on the method
  my $obs;
  if($method =~ /^kf(2)?$/){
    $obs=$refObj->calc_branch_length_distance($queryObj,0);
    # Branch length score is the branch length distance squared
    if(defined($1)){
      $obs=$obs**2;
    }
  } elsif($method eq 'rf'){
    $obs=$refObj->calc_symdiff($queryObj);
    $obs=$obs*2; # this version of RF always returns a single-sided value
  } else {
    die "ERROR: method $method is not implemented";
  }

  die "INTERNAL ERROR: invalid distance score" if(!defined($obs));

  return $obs;
}


sub usage{
  local $0=basename $0;
  "$0: compares two trees
  Usage: $0 ref.dnd query.dnd [query2.dnd...]
  --method    kf     kf:  kuhner-felsenstein (branch length distance)
                     kf2: branch length score (k-f squared)
                     rf:  robinson-foulds (symmetric distance)
                     WARNING: the kf and kf2 metrics do not produce
                     the same values as R-phangorn or Phylip-treedist.
  --numtrees  1000   How many random trees to compare against?
                     Use 0 to not run a statistical test.
  "
}

