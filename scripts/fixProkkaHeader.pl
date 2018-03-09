#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Perl;
use File::Basename qw/fileparse basename dirname/;

exit main();

sub main{
  die usage() if(!@ARGV);
  
  my($infile)=@ARGV;
  my $in=Bio::SeqIO->new(-file=>$infile,-verbose=>-1);
  my $out=Bio::SeqIO->new(-format=>"genbank");
  my $i=0;
  while(my $seq=$in->next_seq){
    my $id=sprintf("contig%06d",++$i);
    $seq->id($id);
    $out->write_seq($seq);
  }
  return 0;
}

sub usage{
  local $0=basename $0;
  "$0: Fixes headers in genbank files
  Usage: $0 in.gbk > out.gbk
  "
}
