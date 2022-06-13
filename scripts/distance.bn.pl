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

# Make a global hash with deflines
my %globalDefline;
# Keep track of the order in which they were added
my @deflineAddedOrder;
# Only one thread at a time can add/delete from this cache using this shared var
my $deflineLock :shared;

local $0 = basename $0;
sub logmsg{local $|=1; my $tid=threads->tid; local $0=basename $0; print STDERR "$0 (TID:$tid): @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help alleles=i key=i mem=i metric=s numcpus=i debug)) or die $!;
  usage() if($$settings{help} || @ARGV < 1);
  $$settings{numcpus} ||= 1;
  $$settings{mem} ||= 0;
  $$settings{key} ||= 1;
  $$settings{alleles} ||= 2;

  if($$settings{metric}){
    logmsg "WARNING: --metric is not configured in this script";
    ...;
  }
  if($$settings{numcpus} > 1){
    logmsg "WARNING: --numcpus is not used in this script";
  }

  my @tsv = @ARGV;
  my %allele;
  for my $file(@tsv){
    my $suballeles = readAlleleTsv($file, $settings);
    # merge allele hashes
    %allele = (%allele, %$suballeles);
  }
  
  # Create a set of comparisons
  logmsg "Setting up the comparisons";
  my @comparison;
  my @fasta = sort {$a cmp $b} keys(%allele);
  my $numFasta = @fasta;
  for(my $i=0; $i<$numFasta; $i++){
    for(my $j=$i+1; $j<$numFasta; $j++){
      push(@comparison,
        [sort {$a cmp $b} ($fasta[$i], $fasta[$j]) ]
      );
    }
  }
  my $numComparisons = scalar(@comparison);
  logmsg "$numComparisons comparisons set up.";

  # Start off printing the table
  print join("\t", qw(sample1 sample2 identity numSame numCompared))."\n";
  my $i=0;
  for my $c(@comparison){
    # Get the distances for this comparison
    my $dist = allelicDistance($allele{$$c[0]}, $allele{$$c[1]}, $settings);
    print join("\t", $$c[0], $$c[1], $$dist{identity}, $$dist{numSame}, $$dist{numCompared});
    print "\n";

    if(++$i % 10000 == 0){
      logmsg "$i distances calculated";
    }
  }
  
  return 0;
}

sub allelicDistance{
  my($alleles1, $alleles2, $settings) = @_;
  
  my @locus = keys(%$alleles1);
  my $numLoci=@locus;
  my $numSame = 0;
  my $numCompared = 0;
  for(my $i=0; $i<$numLoci; $i++){
    # Will be using the value of the alleles multiple times
    # and so might as well pull them out into $a1, $a2
    my $a1 = $$alleles1{$locus[$i]};
    my $a2 = $$alleles2{$locus[$i]};

    # Many loci have been translated to -1 in readAlleleTsv()
    # and if they are -1, we do not want to compare.
    if($a1 == -1 || $a2 == -1){
      next;
    }

    # If we are comparing, tally it
    $numCompared++;
    # If they are the same alleles, tally it
    if($a1 == $a2){
      $numSame++;
    }
  }

  my $dist = {
    numCompared => $numCompared,
    numSame     => $numSame,
    identity    => $numSame/$numCompared,
  };
  return $dist;
}


sub readAlleleTsv{
  my($tsv, $settings) = @_;

  my %allele;
  open(my $fh, $tsv) or die "ERROR: could not read tsv $tsv: $!";
  my $header = <$fh>;
  chomp($header);
  my (@header) = split(/\t/, $header);
  my @alleleName = splice(@header, $$settings{alleles}-1);
  #splice(@header, 0, $$settings{alleles}-1); # remove $key .. $corePercentCalled
  while(my $line = <$fh>){
    chomp $line;
    # Convert ? cells to -1 for purposes of comparison
    $line =~ s/\?/-1/g;

    my @F = split(/\t/, $line);
    my $key = $F[$$settings{key}-1];

    my %F;
    my @alleleInt = splice(@F, $$settings{alleles}-1);
    @F{@alleleName} = @alleleInt;

    $allele{$key} = \%F;

    #if(keys(%allele) > 5){
    #  logmsg "DEBUG";
    #  last;
    #}

    if(keys(%allele) % 500 == 0){
      logmsg "Read ".keys(%allele)." entries";
      #logmsg "DEBUG"; last;
    }
  }

  return \%allele;
}

sub usage{
  print "$0: pairwise distance between two MLST results from BioNumerics 
  Usage: $0 [options] */results_alleles.tsv
  OPTIONS
  --metric  which distance metric to use (default: identity)
  --debug   Print useful debugging messages
  --numcpus Default:1. Tentative sweet spot: 8
  --key     1-based field int where the unique sample name is
  --alleles 1-based field int where alleles start
  --help    This useful help menu
";
  exit 0;
}
