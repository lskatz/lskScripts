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
use List::Util qw/shuffle sum/;

use Bio::Phylo::IO;
use Bio::Phylo::Factory;

use threads;
use Thread::Queue;

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
  GetOptions($settings,qw(help numobserved=i observed-dir=s numtaxa=i background-file=s tempdir=s method=s numcpus=i numtrees=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{numtrees}||=0;
  $$settings{method}||="kf";
  $$settings{method}||=lc($$settings{method});
  $$settings{tempdir}||=tempdir(basename($0).".XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{tempdir}=File::Spec->rel2abs($$settings{tempdir});
  $$settings{numtaxa}||=0;
  $$settings{numobserved}||=1;

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

  # Do the trees have the same taxa, etc?
  validateTrees($ref,\@query,$settings);

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

sub validateTrees{
  my($ref,$queryList,$settings)=@_;

  logmsg "Validating trees...";

  my $errorCount=0;

  my $refObj=Bio::Phylo::IO->parse(
    -format=>"newick",
    -file=>$ref,
  )->first;

  my $taxa=$refObj->get_terminals();
  my %taxonId;
  for(@$taxa){
    $taxonId{$_->get_name}=1;
  }

  for my $query(@$queryList){
    my $queryObj=Bio::Phylo::IO->parse(
      -format=>"newick",
      -file=>$query,
    )->first;

    my @queryTaxonIds=map{$_->get_name} @{ $queryObj->get_terminals() };

    # Check that all query IDs are in ref
    my %queryTaxonId; # make an index
    for my $queryTaxonId(@queryTaxonIds){
      $queryTaxonId{$queryTaxonId}=1; # make an index
      if(!$taxonId{$queryTaxonId}){
        logmsg "ERROR: taxon $queryTaxonId is not in ref tree $ref";
        $errorCount++;
      }
    }
    # Check that all ref IDs are in query
    for my $refTaxonId(keys(%taxonId)){
      if(!$queryTaxonId{$refTaxonId}){
        logmsg "ERROR: taxon $refTaxonId is not in query tree $query";
        $errorCount++;
      }
    }
  }

  die "Exiting due to $errorCount errors" if($errorCount);

  logmsg "Trees have been validated";

  return 1;
}

sub bioPhyloDist{
  my($query,$ref,$method,$settings)=@_;

  my $observedQ = Thread::Queue->new;
  my $printTreeQueue = threads->new(sub{
    my $outdir = $$settings{'observed-dir'};
    if(!defined($outdir)){
      return;
    }

    mkdir $outdir;
    my $outfile = $outdir."/".basename($query);
    open(my $fh, ">", $outfile) or die "ERROR writing to $outfile: $!";
    while(defined(my $tree = $observedQ->dequeue)){
      print $fh $tree."\n";
    }
    close $fh;
  });

  # First find the observed value.
  # Multithreaded so that we can get multiple observed values.
  my @thr;
  for(0..$$settings{numcpus}-1){
    my %settingsCopy=%$settings;
    my $repsPerThread = int($$settings{numobserved}/$$settings{numcpus}) + 1;
    $thr[$_] = threads->new(sub{
      my($repsPerThread, $ref, $query, $observedQ, $settings)=@_;
      my @obs;
      for my $obsRep(1..$repsPerThread){
        my $obs;
        my $tries=0;
        if($ref eq $query){
          push(@obs, 0);
          next;
        }
        
        do{
          my($tmp,$queryDnd,$refDnd)=observedScore($query,$ref,$method,$settings);
          $obs = $tmp;
          $observedQ->enqueue($queryDnd);
          if($tries++ > 10000){
            die "ERROR: tried to get a zero self vs self score for $ref but couldn't after $tries tries. The tree might be too multifurcating. Force bifurcation and try again";
          }
        } while($ref eq $query && $obs != 0);
        push(@obs, $obs);
        logmsg "I calculated the distances between trees $obsRep times";
      }
      return \@obs;
    },$repsPerThread, $ref, $query, $observedQ, \%settingsCopy);
  }

  my @obs;
  for(@thr){
    my $obsTmp = $_->join;
    push(@obs, @$obsTmp);
  }
  $observedQ->enqueue(undef); # kill the printer thread
  @obs = sort {$a<=>$b} (shuffle(@obs))[0..$$settings{numobserved}-1];
  my $obs = sprintf("%0.2f",sum(@obs)/scalar(@obs));

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
    if($$settings{'background-file'}){
      logmsg "Writing raw scores to $$settings{'background-file'}";
      open(my $fh,">", $$settings{'background-file'}) or die "ERROR: could not write to ".$$settings{'background-file'}.": $!";
      print $fh join("\n",sort{$a <=> $b} @scores);
      print $fh "\n";
      close $fh;
    }

    my $stat=Statistics::Descriptive::Full->new();
    $stat->add_data(@scores);
    $avg=$stat->mean;
    $stdev=$stat->standard_deviation;
    $stdev=1e-9 if($stdev <= 0);
    $Z=($obs - $avg)/$stdev;
    $p=cdf($obs,$avg,$stdev);
  }
  
  $printTreeQueue->join;
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
    my $dist=observedScore($randQueryTree,$ref,$method,$settings);
    push(@dist,$dist);
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

  # Choose some random taxa to remove/keep
  if(my $numTaxa = $$settings{numtaxa}){
    my $numLeaves = $queryObj->get_ntax();
    my $numToRemove = $numLeaves - $$settings{numtaxa};
    my @leaf_removal = (
      map { $_->get_name() }
      shuffle(
      @{ $queryObj->get_terminals }
    ))[0..$numToRemove -1];
    $queryObj = $queryObj->prune_tips(\@leaf_removal);
    $refObj   = $refObj->prune_tips(\@leaf_removal);
  }

  # Tree operations for both
  for($queryObj,$refObj){
    # Make polytomies out of short branch lengths
    $_->visit_depth_first(
      -pre=>sub{
        my ($node) = @_;
        my $branch_length = $node->get_branch_length || 0;
        if($branch_length < 1e-5 
          && defined(my $parent=$node->get_parent) 
          && defined(my $grandparent = $node->get_parent->get_parent)
        ){
          my $parentBranchLength = $parent->get_branch_length || 0;
          $branch_length += $parentBranchLength;
          $node->set_branch_length($branch_length);
          $node->set_parent($grandparent);
          #logmsg "Set the parent ". $node->get_name." with length ".$node->get_branch_length;
        }
      },
    );

    # Break polytomies and make the tree binary.
    $_=$_->resolve()->deroot();
  }

  # Set the tree objects equal to each other if the 
  # file paths are the same. Therefore the resolve()
  # function would result in the same tree.
  if($query eq $ref){
    $queryObj=$refObj;
  }

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
  } elsif($method eq 'nrf'){
    $obs=$refObj->calc_symdiff($queryObj);
    $obs=$obs*2; # this version of RF always returns a single-sided value
    $obs=$obs / (2*($refObj->get_ntax() - 3));
  } else {
    die "ERROR: method $method is not implemented";
  }

  die "INTERNAL ERROR: invalid distance score" if(!defined($obs));

  if(wantarray){
    return($obs, $queryObj->to_newick, $refObj->to_newick);
  }
  return $obs;
}


sub usage{
  local $0=basename $0;
  "$0: compares two trees
  Usage: $0 ref.dnd query.dnd [query2.dnd...]
  --method    kf     kf:  kuhner-felsenstein (branch length distance)
                     kf2: branch length score (k-f squared)
                     rf:  robinson-foulds (symmetric distance)
                     nrf: normalized robinson-foulds (symmetric distance)
                     WARNING: the kf and kf2 metrics do not produce
                     the same values as R-phangorn or Phylip-treedist.
  --numobserved  1   Repeat the calculation between the observed
                     trees (not the random trees) and report the average.
  --numtrees  1000   How many random trees to compare against?
                     Use 0 to not run a statistical test.
  --numcpus      1
  --observed-dir     A directory for recording the observed trees
  --background-file  A filename for recording the background
                     distribution distances (optional)
  --numtaxa      0   If > 0, only keep X leaves on the tree for
                     comparisons.
  "
}

