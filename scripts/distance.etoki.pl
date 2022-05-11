#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use List::MoreUtils qw/uniq/;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help metric=s debug)) or die $!;
  usage() if($$settings{help} || @ARGV < 2);

  if($$settings{metric}){
    logmsg "WARNING: --metric is not configured in this script";
  }

  my @fasta = @ARGV;
  print join("\t", qw(fasta1 fasta2 dist))."\n";
  for(my $i=0;$i<@fasta;$i++){
    for(my $j=$i+1; $j<@fasta; $j++){
      my $dist = percentIdentical($fasta[$i], $fasta[$j], $settings);
      print join("\t", 
        # sort genome1 and genome2 to make the output more predictable
        sort({$a cmp $b } ($fasta[$i], $fasta[$j])) ,
        $dist
      );
      print "\n";
    }
  }

  return 0;
}

sub percentIdentical{
  my($fasta1, $fasta2, $settings) = @_;

  my $deflines1 = readFastaDeflines($fasta1, $settings);
  my $deflines2 = readFastaDeflines($fasta2, $settings);

  my $numLociInCommon=0;
  my $numAllelesInCommon=0;
  my @allLoci = sort {$a cmp $b} uniq(keys(%$deflines1), keys(%$deflines2));
  for my $locus(@allLoci){
    if(!defined($$deflines1{$locus}) || !defined($$deflines2{$locus})){
      logmsg "NOTE: Skipping locus $locus" if($$settings{debug});
      next;
    }
    $numLociInCommon++;
    if($$deflines1{$locus}{id} eq $$deflines2{$locus}{id}){
      $numAllelesInCommon++;
    }
  }

  my $percentIdentical = $numAllelesInCommon / $numLociInCommon;
  return $percentIdentical;
}

sub readFastaDeflines{
  my($fasta, $settings) = @_;

  my %defline;

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

  return \%defline;
}

sub usage{
  print "$0: pairwise distance between two MLST results from EToKi
  Usage: $0 [options] *.etoki.fasta
    where etoki.fasta files are EToKi results with MLST deflines
  --metric which distance metric to use (default: identity)
  --debug  Print useful debugging messages
  --help   This useful help menu
";
  exit 0;
}
