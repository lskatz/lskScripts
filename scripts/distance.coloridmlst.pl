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
use Storable qw/dclone/;

# Make a global hash with deflines
my %globalDefline;
# Keep track of the order in which they were added
my @deflineAddedOrder;
# Only one thread at a time can add/delete from this cache using this shared var
my $printLock :shared;

local $0 = basename $0;
sub logmsg{local $|=1; my $tid=threads->tid; local $0=basename $0; print STDERR "$0 (TID:$tid): @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help metric=s numcpus=i debug)) or die $!;
  usage() if($$settings{help} || @ARGV < 1);
  $$settings{numcpus} ||= 1;

  if($$settings{metric}){
    logmsg "WARNING: --metric is not configured in this script";
    ...;
  }

  # Get all allele calls into memory
  my @tsv = @ARGV;
  my %allele;
  for my $file(@tsv){
    my $suballeles = readAlleleTsv($file, $settings);
    # merge allele hashes
    %allele = (%allele, %$suballeles);
  }
  
  # Set up the threads
  my $Q = Thread::Queue->new;
  my @thr;
  for(my $i=0;$i<$$settings{numcpus};$i++){
    $thr[$i] = threads->new(\&allelicDistanceWorker, \%allele, $Q, $settings);
  }
  # Start off printing the table
  print join("\t", qw(sample1 sample2 identity numSame numCompared))."\n";

  # Create a set of comparisons
  # This should be a very fast step and so do not enqueue until the end,
  # since enqueue() is a slow step.
  logmsg "Setting up the comparisons";
  my @comparison;
  my @fasta = sort {$a cmp $b} keys(%allele);
  my $numFasta = @fasta;
  for(my $i=0; $i<$numFasta; $i++){
    for(my $j=$i+1; $j<$numFasta; $j++){
      # Each element of the @comparison is an array ref where
      # the first element comes first alphabetically.
      push(@comparison,
        [sort {$a cmp $b} ($fasta[$i], $fasta[$j]) ]
      );
    }
  }
  my $numComparisons = scalar(@comparison);
  logmsg "$numComparisons comparisons set up.";
  $Q->enqueue(@comparison);
  # Since there will be no more comparisons generated after this,
  # also send the termination signal to the threads.
  # Send a multiplier of how many are batch-dequeued times the number of threads.
  my @undef = (undef) x scalar(@thr);
  $Q->enqueue(@undef);
  #logmsg "Enqueued termination signal by sending ".scalar(@undef)." terms";

  # Wait for threads to finish
  for(@thr){
    # logmsg "Joining TID".$_->tid;
    $_->join;
  }

  return 0;
}

sub allelicDistanceWorker{
  my($allelesOrig, $Q, $settings) = @_;

  my $alleles = dclone($allelesOrig);

  # Counter for status updates
  my $numCompared = 0;

  # locking and printing is a slow step, send 10000 lines for 
  # printing at a time.
  my @distBuffer = ();

  while(my $c = $Q->dequeue){
    my $dist = allelicDistance($$alleles{$$c[0]}, $$alleles{$$c[1]}, $settings);
    push(@distBuffer, 
      join("\t", $$c[0], $$c[1], $$dist{identity}, $$dist{numSame}, $$dist{numCompared})
    );
    next if(@distBuffer < (10 ** 7));

    lock($printLock);
    $numCompared += scalar(@distBuffer);
    logmsg "Sent $numCompared comparisons to print from this thread";
    for(@distBuffer){
      print $_ ."\n";
    }
    @distBuffer = ();
  }
  lock($printLock);
  for(@distBuffer){
    print $_ ."\n";
  }
  $numCompared += scalar(@distBuffer);
  @distBuffer = ();
  logmsg "DONE! Sent $numCompared comparisons to print from this thread";
}

sub allelicDistance{
  my($alleles1, $alleles2, $settings) = @_;
  
  my @locus = uniq(keys(%$alleles1), keys(%$alleles2));
  my $numLoci=@locus;
  my $numSame = 0;
  my $numCompared = 0;
  for(my $i=0; $i<$numLoci; $i++){
    # Will be using the value of the alleles multiple times
    # and so might as well pull them out into $a1, $a2
    # If they are not defined, then we will have 'or next'
    my $a1 = $$alleles1{$locus[$i]} or next;
    my $a2 = $$alleles2{$locus[$i]} or next;

    # If we are comparing, tally it
    $numCompared++;
    # If they are the same alleles, tally it
    if($a1 == $a2){
      $numSame++;
    }
  }

  my $identity = 0;
  if($numCompared > 0){
    $identity = $numSame/$numCompared;
  }

  my $dist = {
    numCompared => $numCompared,
    numSame     => $numSame,
    identity    => $identity,
  };
  return $dist;
}


sub readAlleleTsv{
  my($tsv, $settings) = @_;

  my %allele;
  open(my $fh, $tsv) or die "ERROR: could not read tsv $tsv: $!";
  while(my $line = <$fh>){
    chomp $line;

    my ($hit, $sample, $kmerLength, $percentCoverage) = split(/\t/, $line);
    # Do not allow novel alleles, marked by *
    next if($hit =~ /\*/);

    # Parse locus, allele from hit
    my $lastFieldIdx = rindex($hit, "_");
    my $locus = substr($hit, 0, $lastFieldIdx);
    my $allele = substr($hit, $lastFieldIdx+1);

    # If the locus has been found already, it is a duplicate and so
    # mark it as absent for the purposes of this comparison
    if(defined($allele{$sample}{$locus})){
      $allele{$sample}{$locus} = 0;
    } else {
      $allele{$sample}{$locus} = $allele;
    }

    if($. % 100000 == 0){
      logmsg "Read $. sample/locus lines";
    }
  }
  logmsg "Done reading $. sample/locus lines";
  close $fh;

  return \%allele;
}

sub usage{
  print "$0: pairwise distance between MLST results from ColorID, using alleles.tsv
  Usage: $0 [options] alleles.tsv [alleles2.tsv...]
  OPTIONS
  --metric  which distance metric to use (default: identity)
  --debug   Print useful debugging messages
  --numcpus Default:1. Tentative sweet spot: 8
  --help    This useful help menu
";
  exit 0;
}
