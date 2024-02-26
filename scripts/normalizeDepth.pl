#!/usr/bin/env perl 

use warnings;
use strict;a
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use List::Util qw/shuffle min sum/;
use POSIX qw/ceil/;

use version 0.77;
our $VERSION = '0.2.0';

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help iterations|iters=i paired-end k=i target-depth=i)) or die $!;
  usage() if($$settings{help});
  my $k = $$settings{k} || 11;
  my $pairedEnd = $$settings{'paired-end'} || 0;
  my $iterations = $$settings{iterations} || 1;
  my $maxDepth = $$settings{'target-depth'} or die "ERROR: need --target-depth";
  if($pairedEnd){
    $maxDepth = ceil($maxDepth / 2);
  }

  my $entries = readFastq($pairedEnd, $settings);
  for my $i(1..$iterations){
    logmsg "Starting round $i with ".scalar(@$entries)." spots";
    $entries = normalize($entries, $maxDepth, $k, $pairedEnd, $i, $settings);
  }

  logmsg "DONE! Finished with ".scalar(@$entries)." spots";
  print $_ for(@$entries);

  return 0;
}

sub normalize{
  my($entries, $maxDepth, $k, $pairedEnd, $iteration, $settings) = @_;

  # kmers => {seq1=>1, {seq2=>1}...}
  my %kmers;
  # sequences => [kmer1, kmer2...]
  my %seqs;

  # Read the input fastq
  my $numEntries = @$entries;
  #my %entries;
  #for my $e(@$entries){
  #  $entries{$e} = 1;
  #}

  for(my $i=0; $i<$numEntries; $i++){
    my $entry = $$entries[$i];
    my (undef, $seq, undef, undef, undef, $seq2) = split(/\n/, $entry);

    # How many kmers to get? Length of the sequence minus k, plus 1
    my $numKmers = length($seq)-$k+1;
    my @kmers_in_seq;
    for(my $i=0; $i<$numKmers; $i++){
      # Get the kmer at this position
      my $kmer = substr($seq, $i, $k);

      # Initialize the kmer list of sequences
      $kmers{$kmer} //= {entries=>{}, count=>0};

      # Save this sequence under this kmer
      $kmers{$kmer}{entries}{$entry} = 1;
      $kmers{$kmer}{count}++;

      push(@kmers_in_seq, $kmer);
    }
    # Save the kmers for the sequence.
    $seqs{$entry} = \@kmers_in_seq;
  }
  undef($entries);

  #if($iteration > 10){
  #  print Dumper $kmers{ATTTACAACAT};
  #  die;
  #}

  # Sort kmers by count,
  # then by alphabetical to keep it stable sort.
  my @sortedKmer = sort{
      $kmers{$b}{count} <=> $kmers{$a}{count}
      || $b cmp $a
    } keys(%kmers);

  my @totalCoverage;
  for my $kmer(@sortedKmer){
    push(@totalCoverage, $kmers{$kmer}{count});
  }
  logmsg "Average coverage ".(sum(@totalCoverage)/@totalCoverage);
  logmsg "Kmers with most coverage";
  for(my $j=0;$j<5;$j++){
    logmsg "  $sortedKmer[$j] ".$kmers{$sortedKmer[$j]}{count};
  }
  for(my $j=-5;$j<0;$j++){
    logmsg "  $sortedKmer[$j] ".$kmers{$sortedKmer[$j]}{count};
  }

  # which entries to keep after normalizing
  my @keepEntries = ();

  my $numKmersToView = scalar(@sortedKmer); # min(100000, scalar(@sortedKmer));
  for(my $j=0;$j<$numKmersToView;$j++){
    my $kmer = $sortedKmer[$j];

    my $seqs = $kmers{$kmer}{entries};
    my @shuffled = shuffle(keys(%$seqs));
    my $numReads = min($maxDepth, scalar(@shuffled));
    if($kmer eq 'ATTTACAACAT'){
      logmsg "Will reduce $kmer to $numReads entries";
      logmsg "  Found ".scalar(@shuffled)." reads with $kmer";
    }
    for(my $i=0;$i<$numReads;$i++){
      push(@keepEntries, $shuffled[$i]);
      if($kmer eq 'ATTTACAACAT'){
        my($id) = split(/\n/, $shuffled[$i]);
        logmsg "    $i) Added $id";
      }

      # Remove this from the kmer hash now that it has been used
      for my $kmer(@{ $seqs{$shuffled[$i]} }){
        delete($kmers{$kmer}{entries}{$shuffled[$i]});
      }
      #delete($entries{$shuffled[$i]});
    }
  }

  # Add back on whatever is left that wasn't downsampled
  #push(@keepEntries, keys(%entries));

  return [shuffle @keepEntries];
}

# Read fastq from stdin
sub readFastq{
  my($pairedEnd, $settings) = @_;

  my $fh = \*STDIN;
  #logmsg "DEBUG: different input file"; open(my $fh,"zcat tests/unittests/input/SRR27366697.10x.fastq.gz | ") or die "ERROR: could not open test file: $!";

  my @entries;
  while(my $id = <$fh>){
    # Get the rest of the entry
    my($seq, $plus, $qual) = (scalar(<$fh>), scalar(<$fh>), scalar(<$fh>));
    my $entry = "$id$seq$plus$qual";
    if($pairedEnd){
      my($id2, $seq2, $plus2, $qual2) = (scalar(<$fh>), scalar(<$fh>), scalar(<$fh>), scalar(<$fh>));
      $entry .= "$id2$seq2$plus2$qual2";
    }
    push(@entries, $entry);
  }
  close $fh;
  return \@entries;
}


sub usage{
  print "$0: normalizes depth by kmer coverage
  Usage: zcat in.fastq.gz | $0 [options] | gzip -c > out.fastq.gz
  --target-depth  x  The maximum kmer coverage
  --iterations    1  The number of times to downsample
  --paired-end       If interleaved paired end
  --k             11 Kmer length
  --help   This useful help menu
  \n";
  exit 0;
}
