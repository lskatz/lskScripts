#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>
# Splits up a bionumerics fasta file into contigs

use strict;
use warnings;
use Bio::Perl;
use File::Basename;
use autodie;

exit main();

sub main{
  die usage() if(!@ARGV);
  
  my $in=Bio::SeqIO->new(-file=>$ARGV[0]); 
  my $out=Bio::SeqIO->new(-format=>"fasta"); 
  while(my $seq=$in->next_seq){
    my @seq=split(/\|/,$seq->seq);
    my $id=$seq->id;
    for(my $i=1;$i<=@seq;$i++){
      my $subseq=Bio::Seq->new(-seq=>$seq[$i-1],-id=>$id."_".$i);
      $out->write_seq($subseq);
    }
  }

  return 0;
}

sub usage{
  local $0=fileparse $0;
  "Usage: $0 assembly.fasta > out.fasta
  "
}
