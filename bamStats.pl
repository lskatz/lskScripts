#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  die usage() if($$settings{help});
  die usage() if(-t STDIN);

  my $QC=bamMetrics($settings);

  my @header=sort {$a cmp $b} keys(%$QC);

  print join("\t",@header)."\n";
  for my $header(@header){
    print $$QC{$header}."\t";
  }
  print "\n";

  return 0;
}

sub bamMetrics{
  my($settings)=@_;

  my %QC;

  while(<>){
    chomp;
    my($seqid, $flag, $rname, $pos, $mapQ, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = split(/\t/, $_);

    # Individual tags
    if($flag & 4){
      $QC{unmapped}++;
    }
    if(! $flag & 2){
      $QC{improperPair}++;
    }
    
    # Combination tags
    if($flag =~ /^(?:73|133|89|121|165|181|101|117|153|185|69|137)$/){
      $QC{singletonMap}++;
    } elsif ($flag =~ /^(?:77|141)$/){
      $QC{bothUnmapped}++;
    } elsif($flag =~ /^(?:99|147|83|163)$/){
      $QC{bothProperPair}++;
    } elsif($flag =~ /^(?:67|131|115|179)$/){
      $QC{wrongOrientation}++;
    } elsif($flag =~ /^(?:81|161|97|145|65|129|113|177)$/){
      $QC{wrongInsertSize}++;
    }
  }
  return \%QC;
}

sub usage{
  "$0: get QC information on a sam file
  Usage: samtools view file.bam | $0
  "
}
