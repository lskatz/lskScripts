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
  GetOptions($settings,qw(help database=s)) or die $!;
  usage() if($$settings{help} || @ARGV < 1);
  $$settings{mem} ||= 0;
  $$settings{database} || die "ERROR: need --database";

  if($$settings{metric}){
    logmsg "WARNING: --metric is not configured in this script";
    ...;
  }

  my $loci = lociFromDatabase($$settings{database}, $settings);

  # Print off the header with assembly name and loci
  print join("\t", "assembly", @$loci)."\n";

  # Print each row of the matrix, one per assembly
  my @fasta = @ARGV;
  for(my $i=0; $i<@fasta; $i++){
    my $gappedMatrix = gappedMatrix($fasta[$i], $loci, $settings);
    print $gappedMatrix . "\n";
  }

  return 0;
}

sub lociFromDatabase{
  my($db, $settings) = @_;

  my %locus;

  open(my $fh, $db) or die "ERROR: could not read database $db: $!";
  while(<$fh>){
    chomp;
    my ($md5, $locus, $allele) = split /,/;
    $locus{$locus}++;
  }
  close $fh;

  my @locus = sort{$a cmp $b} keys(%locus);
  return \@locus;
}

sub gappedMatrix{
  my($fasta, $loci, $settings) = @_;

  # Initialize the matrix row with the name of the fasta file
  my $row = "$fasta";

  my $deflines = readFastaDeflines($fasta, $settings);

  # Loci are already sorted and so do not sort here
  my $numLoci = @$loci;
  for(my $i=0; $i<$numLoci; $i++){
    my $id = $$deflines{$$loci[$i]}{id};
    if(!defined($id)){
      $id = -1;
    }
    #if($id ne "-1" && $id =~ /\D/){
    #  logmsg "id was $id and now it is -1";
    #  $id = -1;
    #}

    $row .= "\t".$id;
  }
  return $row;
}


sub readFastaDeflines{
  my($fasta, $settings) = @_;

  my %defline;

  # Check if the defline is in the cache and if so,
  # grab it. Lock this step just in case the cache
  # is being modified elsewhere.
  {
    lock($deflineLock);
    if(defined $globalDefline{$fasta}){
      return $globalDefline{$fasta};
    }
  }

  open(my $fh, $fasta) or die "ERROR reading $fasta: $!";
  while(<$fh>){
    # only read deflines for this operation
    next if(!/^>/);
    # Remove whitespace and the initial >
    s/^\s+|\s+$|^>//g;

    my %F;
    my @F = split(/\s+/, $_);
    my $id = shift(@F);
    for my $keyvalue(@F){
      my($key, $value) = split(/=/, $keyvalue, 2);
      $F{$key} = $value;
    }
    # Some notes on the accepted field
    # 1: accepted
    # 2: low-quality
    # 4, 8, 16: Not used in standalone version
    # 32: multiple copies in the genome
    # 64: fragmented allele
    # 128: overlapped with other loci
    # 256: low identity
    # 34: ????
    # 66: ????
    if($F{accepted} > 2){
      logmsg "For locus $id, the accepted field is $F{accepted} and so I will not save this locus for $fasta" if($$settings{debug});
      next;
    }

    if(defined $defline{$id}){
      logmsg "WARNING: duplicate identifier $id found in $fasta";
    }
    $defline{$id} = \%F;
  }
  close $fh;

  # If %globalDefline is too large, delete an entry
  {
    lock($deflineLock);
    if($$settings{mem} > 0 && scalar(keys(%globalDefline)) > $$settings{mem}){
      my $toDelete = shift(@deflineAddedOrder);
      delete($globalDefline{$toDelete});
      logmsg "Removed from cache: $toDelete" if($$settings{debug});
    }
  }

  # Remember these deflines to avoid costly future file IO costs 
  # Set the lock in case the cache is being acted upon
  {
    lock($deflineLock);
    $globalDefline{$fasta} = \%defline;
    logmsg "Added to cache: $fasta" if($$settings{debug});
    push(@deflineAddedOrder, $fasta);
  }

  return \%defline;
}

sub usage{
  print "$0: Generate a matrix of alleles from a set of etoki fasta results
  Usage: $0 [options] *.etoki.fasta > matrix.tsv
    where etoki.fasta files are EToKi results with MLST deflines
  OPTIONS
  --database  The path to the etoki database whose format is CSV with
              three columns: md5sum, locus, allele [required]
  --help    This useful help menu
";
  exit 0;
}
