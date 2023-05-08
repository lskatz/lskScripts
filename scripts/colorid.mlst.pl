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
  $$settings{tempdir} //= tempdir("$0.XXXXXX", CLEANUP=>1, TMPDIR=>1);
  $$settings{k}       ||= 39;
  $$settings{k} >= 3 || die "ERROR: --k should be >= 3";
  $$settings{numcpus} ||= 1;

  which("colorid") or die "ERROR: could not find colorid in your path!";

  my $mlstDir = shift(@ARGV);

  #logmsg "Reading the MLST directory";
  #my $mlstFasta = readMlstFastas($mlstDir, $settings);
  #die system("head $mlstFasta");

  logmsg "Making samples file for building the bigsi db";
  my $samplesFile = buildSamples(\@ARGV, $settings);
  logmsg "Making the bigsi db";
  my $db = buildDb($samplesFile, $settings);
  logmsg "Searching the MLST database";
  my $search = mlst($mlstDir, $db, $settings);
  logmsg "Transforming hits into alleles";
  my $profiles = profiles($search, $settings);

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


sub mlst{
  my($mlstDir, $bxi, $settings) = @_;

  my $results = "$$settings{tempdir}/hits.tsv";

  my @fasta = glob("$mlstDir/*.fasta");

  # Chunk the fasta files into the queue
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

  logmsg "  Waiting on threads, and then I will merge the results."
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

sub mlstWorker{
  my($Q, $bxi, $settings) = @_;

  my $TID = threads->tid;

  my $results = "$$settings{tempdir}/hits.$TID.tsv";
  my $log     = "$$settings{tempdir}/hits.$TID.log";

  open(my $resultsFh, ">", $results) or die "ERROR: could not write to $results: $!";

  while(defined(my $fastas = $Q->dequeue)){
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

sub profiles{
  my($hits, $settings) = @_;

  my $profilesFile = "$$settings{tempdir}/profiles.tsv";

  my %allele;
  my %locusIndex;
  open(my $fh, $hits) or die "ERROR: could not open $hits: $!";
  while(<$fh>){
    chomp;
    my($locusAllele, $sample, $length, $identity) = split /\t/;

    my($locus, $allele) = $locusAllele =~ /(.+)_(.+)/;
    $allele{$sample}{$locus} = $allele;
    $locusIndex{$locus}++;
  }
  close $fh;

  open(my $outFh, ">", $profilesFile) or die "ERROR: could not write to $profilesFile: $!";

  my @sortedLocus = sort{$a cmp $b} keys(%locusIndex);
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


# Read the cgMLST directory and return a single fasta file
sub readMlstFastas{
  my($mlstDir, $settings)=@_;

  my $newFasta = "$$settings{tempdir}/query.fasta";
  open(my $fh, ">", $newFasta) or die "ERROR: writing to $newFasta: $!";

  for my $f(glob("$mlstDir/*.fasta")){

    # Read this file, sequence by sequence
    my ($n, $slen, $comment, $qlen) = (0, 0, 0);
    my @aux = undef;
    my ($id, $sequence);
    open(my $inFh, $f) or die "ERROR: could not read $f: $!";
    while ( ($id, $sequence, undef) = readfq($inFh, \@aux)) {
      next if(length($sequence) < $$settings{k});
      # Avoid low complexity homopolymers
      next if($sequence =~ /A{30,}|C{30,}|G{30,}|T{30,}/);
      # Avoid lots of ambiguities (5+ Ns)
      next if($sequence =~ /N{5,}/);
      print $fh ">$id\n$sequence\n";
    }
    close $inFh;
  }
  close $fh;

  return $newFasta;
}

# Read fq subroutine from Andrea which was inspired by lh3
sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!(@$aux)); # remove deprecated 'defined'
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my $name = /^.(\S+)/? $1 : '';
    my $comm = /^.\S+\s+(.*)/? $1 : ''; # retain "comment"
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $seq, $comm, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq, $comm);
}

sub usage{
  print "$0: runs colorid/mlst on a sample
  Usage: $0 [options] MLST.db/ *.fasta *.fastq.gz
  MLST.db     The MLST database of fasta files, one fasta per locus
  *.fasta *.fastq.gz can be any number of sequence files

  -k          kmer length [default: 39]
  --tempdir   Alternative location for temporary directory
  --help      This useful help menu
  \n";
  exit 0;
}

