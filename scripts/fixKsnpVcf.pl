#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Perl;
use Getopt::Long;
use Data::Dumper;
use constant reportEvery=>1000;

sub logmsg{print STDERR "@_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(ref|reference=s help)) or die $!;
  die usage() if($$settings{help} || !$$settings{ref});

  my %seq;
  my $in=Bio::SeqIO->new(-file=>$$settings{ref});
  while(my $seq=$in->next_seq){
    $seq{$seq->id}=uc($seq->seq);
  }
  
  my $lineCount=0;
  my $numFixed=0;
  while(<>){
    # Print headers
    if(/^#/){
      print;
      next;
    }

    # Fix VCF lines
    $lineCount++;
    chomp;
    my @F=split /\t/;
    $F[1]||=0;
    $F[2]=uc($F[2]);

    # Fix up the kmer line and print it
    $numFixed += !! fixPosition(\@F,\%seq,$settings);

    if($lineCount % reportEvery == 0){
      my $percent=int($numFixed/$lineCount * 100);
      logmsg "Have fixed $numFixed out of $lineCount ($percent%)";
    }
  }

  my $percent=int($numFixed/$lineCount * 100);
  logmsg "Fixed $numFixed out of $lineCount ($percent%)";
}

sub fixPosition{
  my($F,$seqHash,$settings)=@_;

  # Remove the old Chrom/pos information from the VCF line.
  # It's not what we want anyway.
  splice(@$F,0,2);

  # Figure out where the snp is within the kmer to help
  # with the genomic position later on.
  my $dotIndex=index($$F[0],'.');
  # Also need to search with the reverse complement.
  my $revcom=revcom($$F[0])->seq;

  # Keep track of whether or not there are matches
  my $numMatches=0;
  for my $id(keys(%$seqHash)){
    while($$seqHash{$id}=~/($$F[0]|$revcom)/g){
      my $pos=length($`)+$dotIndex;
      print join("\t",$id,$pos,@$F)."\n";
      $numMatches++;
    }
  }

  logmsg "WARNING: I could not find kmer $$F[0] in $$settings{ref}" if($numMatches < 1);
  return $numMatches;
}

sub usage{
  "Usage: $0 -ref reference.fasta < ksnp.vcf > fixed.vcf
  "
}
