#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Splits up a bionumerics fasta file into contigs

use strict;
use warnings;
use Bio::Perl;
use Bio::Tools::GuessSeqFormat;
use File::Basename;
use autodie;
use Getopt::Long;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help outdir=s)) or die $!;

  die usage() if(!@ARGV || $$settings{help});

  if($$settings{outdir}){
    mkdir $$settings{outdir};
  }

  # Because we are reading bionumerics files and do not
  # trust their extensions, make a better guess at their
  # format using the format guesser.
  my $formatGuesser=Bio::Tools::GuessSeqFormat->new(-file=>$ARGV[0]);
  my $format       =$formatGuesser->guess;

  my $in=Bio::SeqIO->new(-file=>$ARGV[0],-format=>$format); 
  while(my $seq=$in->next_seq){
    my @seq=split(/\|/,$seq->seq);
    my $id=$seq->id;
       $id=~s/^denovo\|//;   # remove 'denovo' since most exports seem to have that
    print STDERR "$id\n";

    my $out=Bio::SeqIO->new(-format=>"fasta"); 
    if($$settings{outdir}){
      # For potential filenames, get a safe name
      my $id_safe=$id;
         $id_safe=~s/[^\w\d]//g;  # remove non words, non letters

      my $outfile="$$settings{outdir}/$id_safe.fasta";
      $out=Bio::SeqIO->new(-format=>"fasta",-file=>">$outfile");
    }
    for(my $i=1;$i<=@seq;$i++){
      my $subseq=Bio::Seq->new(-seq=>$seq[$i-1],-id=>$id."_".$i);
      $out->write_seq($subseq);
    }
    $out->close;
  }

  return 0;
}

sub usage{
  local $0=fileparse $0;
  "Usage: $0 bionumerics.fasta > out.fasta
  --outdir   ''  If given, all genomes will be written here.
                 If blank, output will be sent to stdout.
  "
}
