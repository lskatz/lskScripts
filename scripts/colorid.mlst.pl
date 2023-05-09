#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use File::Which qw/which/;

use threads;
use Thread::Queue;

use version 0.77;
our $VERSION = '0.1.1';

local $0 = basename $0;
sub logmsg{local $0=basename $0; my $TID=threads->tid; print STDERR "$0($TID): @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help k=i numcpus=i tempdir=s)) or die $!;
  usage() if(@ARGV < 2 || $$settings{help});

  # Set some default variables
  $$settings{tempdir} //= tempdir("$0.XXXXXX", CLEANUP=>1, TMPDIR=>1);
  $$settings{k}       ||= 39;
  # sanity check on k
  $$settings{k} >= 3 || die "ERROR: --k should be >= 3";
  $$settings{numcpus} ||= 1;

  logmsg "Temporary directory will be at $$settings{tempdir}";

  which("colorid") or die "ERROR: could not find colorid in your path!";

  my $mlstDir = shift(@ARGV);

  logmsg "Making samples tsv file for building the bigsi db";
  my $samplesFile = buildSamples(\@ARGV, $settings);
  logmsg "Recording all loci from the schema";
  my $sortedLoci = readLoci($mlstDir, $settings);
  logmsg "Making the bigsi db";
  my $db = buildDb($samplesFile, $settings);
  logmsg "Searching the MLST database";
  my $search = mlst($mlstDir, $db, $settings);
  logmsg "Transforming hits into alleles";
  my $profiles = profiles($search, $sortedLoci, $settings);

  # Print the profiles file to stdout
  open(my $fh, $profiles) or die "ERROR: could not open for reading $profiles: $!";
  while(<$fh>){
    print;
  }
  close $fh;

  return 0;
}

sub buildSamples{
  my($seqs, $settings) = @_;

  my $samplesTsv = "$$settings{tempdir}/samples.tsv";

  open(my $fh, ">", $samplesTsv) or die "ERROR: could not write to $samplesTsv: $!";
  for my $s(@$seqs){
    print $fh join("\t", $s, $s)."\n";
  }
  close $fh;

  return $samplesTsv;
}

# Get a sorted list of loci from the MLST database
sub readLoci{
  my($dir, $settings) = @_;

  my %locus;

  for my $fasta(glob("$dir/*.fasta")){
    open(my $fh, $fasta) or die "ERROR: could not read fasta file $fasta: $!";
    while(<$fh>){
      chomp;
      if(/^>(.+)_(.+?)\s*/){
        my $locus = $1;

        # Don't use loci with timestamps
        next if($locus =~ /:/);

        $locus{$locus}++;
      }
    }
    close $fh;
  }

  my @loci = sort{$a cmp $b} keys(%locus);

  return \@loci;
}


# Build the bigsi database
sub buildDb{
  my($samplesTsv, $settings) = @_;

  my $db = "$$settings{tempdir}/bxi";
  my $buildLog = "$$settings{tempdir}/build.log";

  system("colorid build -b $db -s 30000000 -n 2 -k $$settings{k} -t $$settings{numcpus} -r $samplesTsv 2> $buildLog 1>&2");
  if($?){
    logmsg "ERROR with colorid build. Here is the log:\n".`cat $buildLog`;
    die;
  }

  return "$db.bxi";
}


# Run the colorid search function
# in threads
sub mlst{
  my($mlstDir, $bxi, $settings) = @_;

  my $results = "$$settings{tempdir}/hits.tsv";

  my @fasta = glob("$mlstDir/*.fasta");

  # Chunk the fasta files into the queue
  # Currently enqueuing 100 at a time.
  # This will let us load the bigsi db fewer times per thread.
  my $Q = Thread::Queue->new();
  while(@fasta){
    $Q->enqueue(
      [splice(@fasta, 0, 100)]
    );
  }

  my @thr;
  for(my $i=0;$i<$$settings{numcpus};$i++){
    $thr[$i] = threads->new(\&mlstWorker, $Q, $bxi, $settings);
  }

  # Send termination signals to the threads
  for(@thr){
    $Q->enqueue(undef);
  }

  # Merge all the threads results files into one
  logmsg "  Waiting on threads, and then I will merge the results.";
  open(my $resultFh, ">", $results) or die "ERROR: could not write to larger results file: $!";
  for my $t(@thr){
    local $/=undef;
    my $resultsFile = $t->join;
    open(my $fh, $resultsFile) or die "ERROR: could not read thread results file $resultsFile: $!";
    my $content = <$fh>;
    close $fh;

    print $resultFh $content;
  }
  close $resultFh;

  return $results;
}

# This is the threads for `colorid search`
sub mlstWorker{
  my($Q, $bxi, $settings) = @_;

  my $TID = threads->tid;

  my $results = "$$settings{tempdir}/hits.$TID.tsv";
  my $log     = "$$settings{tempdir}/hits.$TID.log";

  open(my $resultsFh, ">", $results) or die "ERROR: could not write to $results: $!";

  while(defined(my $fastas = $Q->dequeue)){
    # Each dequeue yields a chunk of fasta files
    my $fastaStr = join(" ", @$fastas);
    system("colorid search -b $bxi -q $fastaStr -m -s > $results.tmp 2>$log");
    if($?){
      die "ERROR with colorid search with $fastaStr\n  Here is the log:\n".`cat $log`;
    }

    # Open the output file for this thread and this fasta file
    # and do it in this local context so that the file gets 
    # automatically closed at the end and so that we can
    # set the line terminator to undef so that we just slurp
    # it up.
    # Then, append the file to the overall results file.
    {
      local $/=undef;
      open(my $fh, "$results.tmp") or die "ERROR: could not read $results.tmp: $!";
      my $content = <$fh>;
      print $resultsFh $content;
    }

  }
  close $resultsFh;

  return $results;
}

# Read the hits file from colorid and make a profiles file:
# a file with one line per sample and one locus/allele per col
sub profiles{
  my($hits, $loci, $settings) = @_;

  my $profilesFile = "$$settings{tempdir}/profiles.tsv";

  # Make a 2d hash of sample => {locus=>allele}
  my %allele;
  my %locusIndex;
  open(my $fh, $hits) or die "ERROR: could not open $hits: $!";
  while(<$fh>){
    chomp;
    my($locusAllele, $sample, $length, $identity) = split /\t/;

    # Right split
    my($locus, $allele) = $locusAllele =~ /(.+)_(.+)/;
    $allele{$sample}{$locus} = $allele;
    $locusIndex{$locus}++;
  }
  close $fh;

  # Print to the profiles file
  open(my $outFh, ">", $profilesFile) or die "ERROR: could not write to $profilesFile: $!";

  my @sortedLocus = @$loci;
  my @sortedSample= sort{$a cmp $b} keys(%allele);

  print $outFh join("\t", "sample", @sortedLocus)."\n";
  for my $sample(@sortedSample){
    my $line = $sample;
    for my $locus(@sortedLocus){
      $allele{$sample}{$locus} //= "-1";
      $line .= "\t". $allele{$sample}{$locus};
    }
    $line .= "\n";
    print $outFh $line;
  }

  return $profilesFile;
}


sub usage{
  print "$0: runs colorid/mlst on samples and produces a profiles.tsv in stdout
  Usage: $0 [options] MLST.db/ *.fasta *.fastq.gz > profiles.tsv
  MLST.db     The MLST database of fasta files, one fasta per locus
  *.fasta *.fastq.gz can be any number of sequence files

  -k          kmer length [default: 39]
  --tempdir   Alternative location for temporary directory
  --numcpus   Default: 1
  --help      This useful help menu
  \n";
  exit 0;
}

