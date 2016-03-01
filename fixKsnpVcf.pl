#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Perl;
use Getopt::Long;

sub logmsg{print STDERR "@_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(ref|reference=s only-positional help)) or die $!;
  die usage() if($$settings{help} || !$$settings{ref});

  my %seq;
  my $in=Bio::SeqIO->new(-file=>$$settings{ref});
  while(my $seq=$in->next_seq){
    $seq{$seq->id}=$seq->seq;
  }
  
  while(<>){
    chomp;
    my @F=split /\t/;
    $F[1]||='.';

    if(/^#/){
      
    } else {
      my($contig,$pos)=findPosition($F[2],\%seq,$settings);
      if($contig && $pos){
        $F[0]=$contig;
        $F[1]=$pos;
      } elsif($$settings{'only-positional'}){
        next;
      }
    }

    print join("\t",@F)."\n";
  }
}

sub findPosition{
  my($kmerRegex,$seqHash,$settings)=@_;

  my $dotIndex=index($kmerRegex,'.');
  my @ID=keys(%$seqHash);

  for my $id(@ID){
    my $seq=$$seqHash{$id};
    if($seq=~/($kmerRegex)/i){
      my $pos=length($`)+$dotIndex;
      return ($id,$pos);
    }
  }

  logmsg "ERROR: I could not find kmer $kmerRegex in $$settings{ref}";
  return (undef,undef);
}

sub usage{
  "Usage: $0 -ref reference.fasta < file.vcf > fixed.vcf
  --only-positional  Removes any position whose kmer is not found
                     in the reference fasta
  "
}
