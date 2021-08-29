#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  usage() if($$settings{help} || @ARGV < 2);

  my $refFasta = shift(@ARGV);

  for my $bam (@ARGV){
    printFastq($bam, $refFasta, $settings);
  }

  return 0;
}

sub printFastq{
  my($bam, $refFasta, $settings) = @_;

  open(my $fh, "samtools view '$bam' |") or die "ERROR using samtools view on $bam: $!";
  while(my $line = <$fh>){
    chomp($line);
    my($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) = 
        split(/\t/, $line);

    # if the read is not unmapped then replace with reference
    ...;
    if(! $flag & 4){
      my $refHit = referenceHit($rname, $pos, $cigar, $refFasta, $settings);
      $seq = $refHit;
      $qname .= " replaced";
    }

    print "\@$qname\n$seq\n+\n$qual\n";
  }
}

sub referenceHit{
  my($rname, $pos, $cigar, $refFasta, $settings) = @_;

  # Determine length from cigar
  # TODO other operation codes like [NSHP=X]
  my $length = 0;
  while($cigar =~ /(\d+)(\w)/g){
    my $code = $2;
    my $int  = $1;
    if($code eq 'M'){
      $length+=$int;
    } elsif($code eq 'D') {
      $length+=$int;
    } elsif($code eq 'I') {
      $length+=0;
    } else {
      die "ERROR: cigar string has a $code which I do not know how to interpret.  Here is the full cigar string: $cigar";
    }
  }

  my $stopPos = $pos + $length - 1;
  my $refHit = `samtools faidx $refFasta '$rname:$pos-$stopPos' | tail -n +2`;
  die "ERROR running samtools faidx on $refFasta" if $?;
  chomp($refHit);
  $refHit =~ s/\n//g;
  $refHit =~ tr/[a-z]/[A-Z]/;

  return $refHit;
}

sub usage{
  print "$0: print a bam as a fastq file, replacing reads with the reference genome
  Usage: $0 [options] ref.fasta *.bam > out.fastq
  --help   This useful help menu
  ";
  exit 0;
}
