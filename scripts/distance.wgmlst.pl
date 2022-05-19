#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use List::MoreUtils qw/uniq/;
use threads;
use Thread::Queue;
use threads::shared;

# Make a global hash with allele calls 
my %globalCall;
# Keep track of the order in which they were added
my @callAddedOrder;
# Only one thread at a time can add/delete from this cache using this shared var
my $callLock :shared;

local $0 = basename $0;
sub logmsg{local $|=1; my $tid=threads->tid; local $0=basename $0; print STDERR "$0 (TID:$tid): @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help mem=i metric=s numcpus=i debug)) or die $!;
  usage() if($$settings{help} || @ARGV < 2);
  $$settings{numcpus} ||= 1;
  $$settings{mem} ||= 0;

  if($$settings{metric}){
    logmsg "WARNING: --metric is not configured in this script";
    ...;
  }

  my @wgmlst = @ARGV;
  my @comparison;
  for(my $i=0; $i<@wgmlst; $i++){
    for(my $j=$i+1; $j<@wgmlst; $j++){
      push(@comparison, [$wgmlst[$i], $wgmlst[$j]]);
    }
  }

  # Kick off the threads
  my $comparisonQueue = Thread::Queue->new(@comparison);
  my @thr;
  my $printQueue = Thread::Queue->new();
  for(0..$$settings{numcpus}-1){
    push(@thr,
      threads->new(\&percentIdenticalWorker, $comparisonQueue, $printQueue, $settings));
    # Give this thread an eventual termination signal
    $comparisonQueue->enqueue(undef); 
  }
  my $printerThr = threads->new(\&printer, $printQueue, $settings);

  for(@thr){
    $_->join;
  }
  # Terminator signal to the printer now that the threads have joined
  $printQueue->enqueue(undef);
  $printerThr->join;

  return 0;
}

# A single thread to get only one place printing, to 
# avoid any race conditions in stdout
sub printer{
  my($queue, $settings) = @_;
  
  # Start off the header
  print join("\t", qw(sample1 sample2 dist numAllelesSame numLociSame))."\n";

  while(defined(my $dist = $queue->dequeue)){
    print join("\t", 
      # sort genome1 and genome2 to make the output more predictable
      sort({$a cmp $b }
        ($$dist{sample1}, $$dist{sample2})
      ),
      sprintf("%0.2f",$$dist{identity}),
      $$dist{numAllelesSame},
      $$dist{numLociSame},
    );
    print "\n";
  }
}

sub percentIdenticalWorker{
  my($inQueue, $outQueue, $settings) = @_;

  logmsg "LAUNCHED thread";

  my $numComparisons = 0;

  while(defined(my $twoGenomes = $inQueue->dequeue)){
    my ($genome1, $genome2) = @$twoGenomes;
    my $dist = percentIdentical($genome1, $genome2, $settings);
    $outQueue->enqueue($dist);
    $numComparisons++;
  }

  logmsg "FINISHED thread. $numComparisons comparisons calculated";
}

sub percentIdentical{
  my($res1, $res2, $settings) = @_;

  my $calls1 = readWgmlstResults($res1, $settings);
  my $calls2 = readWgmlstResults($res2, $settings);

  my $numLociInCommon=0;
  my $numAllelesInCommon=0;
  my @allLoci = sort {$a cmp $b} uniq(keys(%$calls1), keys(%$calls2));
  for my $locus(@allLoci){
    if(!defined($$calls1{$locus}) || !defined($$calls2{$locus})){
      logmsg "NOTE: Skipping locus $locus" if($$settings{debug});
      next;
    }
    $numLociInCommon++;
    if($$calls1{$locus}{allele} eq $$calls2{$locus}{allele}){
      $numAllelesInCommon++;
    }
  }

  my $percentIdentical = $numAllelesInCommon / $numLociInCommon;
  my $return = {
    sample1 => $res1,
    sample2 => $res2,
    identity => $percentIdentical, 
    numAllelesSame => $numAllelesInCommon,
    numLociSame => $numLociInCommon,
  };
  return $return;
}

sub readWgmlstResults{
  my($dir, $settings) = @_;

  # Store allele calls
  my %call;

  # Check if these calls are in the cache and if so,
  # grab it. Lock this step just in case the cache
  # is being modified elsewhere.
  {
    lock($callLock);
    if(defined $globalCall{$dir}){
      return $globalCall{$dir};
    }
  }

  # The alleles call format is not documented and has the following fields
  # separated by tab.
  my @allelesField = qw(locus seq queryContig start stop identity alleleLength allele);

  open(my $allelesFh, "$dir/alleles") or die "ERROR: could not read $dir/alleles: $!";
  while(<$allelesFh>){
    chomp;
    my %F;
    @F{@allelesField} = split /\t/;

    # Save on memory
    for my $key(qw(seq seq queryContig start stop identity alleleLength)){
      delete($F{$key});
    }

    if(defined $call{$F{locus}}){
      logmsg "WARNING: found locus $F{locus} twice in $dir";
    }
    $call{$F{locus}} = \%F;
  }
  close $allelesFh;

  # If %globalCall is too large, delete an entry
  {
    lock($callLock);
    if($$settings{mem} > 0 && scalar(keys(%globalCall)) > $$settings{mem}){
      my $toDelete = shift(@callAddedOrder);
      delete($globalCall{$toDelete});
      logmsg "Removed from cache: $toDelete" if($$settings{debug});
    }
  }

  # Remember these calls to avoid costly future file IO costs 
  # Set the lock in case the cache is being acted upon
  {
    lock($callLock);
    $globalCall{$dir} = \%call;
    logmsg "Added to cache: $dir" if($$settings{debug});
    push(@callAddedOrder, $dir);
  }

  return \%call;
}

sub usage{
  print "$0: pairwise distance between two MLST results from NCBI's wgmlst
  Usage: $0 [options] *.wgmlst
    where *.wgmlst files are directories with wgmlst results.
    A wgmlst directory contains files: alleles, mappings, wgmlst.log, wgmlst.out
  OPTIONS
  --metric  which distance metric to use (default: identity)
  --debug   Print useful debugging messages
  --numcpus Default:1. Tentative sweet spot: 8
  --mem     How many wgmlst allele results to keep in memory (default: 0 which keeps all)
  --help    This useful help menu
";
  exit 0;
}
