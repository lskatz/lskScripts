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

# Print the fastq file from the bam
sub printFastq{
  my($bam, $refFasta, $settings) = @_;

  open(my $fh, "samtools sort -n '$bam' | samtools view -f 1 |") or die "ERROR using samtools view on $bam: $!";
  while(my $line = <$fh>){
    chomp($line);
    my($qname, $flag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual) =
        split(/\t/, $line);

    # add /1 or /2
    if($flag & 0x40){
      $qname .= "/1";
    }
    if($flag & 0x80){
      $qname .= "/2";
    }

    # if the read is not unmapped then replace with reference
    if(! ($flag & 0x4)){
      my $refHit = referenceHit($rname, $pos, $cigar, $refFasta, $settings);
      $seq = $refHit;
      $qname .= " replaced";
    }

    # Sanity check to match the lengths of seq and qual
    if(length($seq) != length($qual)){
      my$bp = length($seq) - length($qual);
      logmsg "WARNING: seq is not the same length as qual for $qname($bp bp)";
      # Adjust to match shortest length
      if(length($seq) > length($qual)){
        $seq = substr($seq, 0, length($qual));
      }else{
        $qual = substr($qual, 0, length($seq));
      }

    }

    print "\@$qname\n$seq\n+\n$qual\n";
  }
}

# Get the sequence of the reference genome at the mapped position
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

  if($length < 1){
    die "INTERNAL ERROR: length of reference hit for this mapped read is <1" . Dumper \@_;
  }

  # Grab the reference hit
  my $stopPos = $pos + $length - 1;
  my $refHit = `samtools faidx $refFasta '$rname:$pos-$stopPos' | tail -n +2`;
  die "ERROR running samtools faidx on $refFasta" if $?;
  chomp($refHit);
  $refHit =~ s/\n//g;
  $refHit =~ tr/[a-z]/[A-Z]/; # uppercase

  return $refHit;
}

sub usage{
  print "$0: print a bam as a fastq file, replacing reads with the reference genome
  Usage: $0 [options] ref.fasta *.bam > out.fastq
  --help   This useful help menu
  ";
  exit 0;
}
