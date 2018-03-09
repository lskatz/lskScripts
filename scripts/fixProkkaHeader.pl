#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Perl;
use File::Basename qw/fileparse basename dirname/;
use Getopt::Long;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help min_length|min-length=i)) or die $!;
  die usage() if(!@ARGV || $$settings{help});
  $$settings{min_length}||=1;
  
  my($infile)=@ARGV;
  my $in=Bio::SeqIO->new(-file=>$infile,-verbose=>-1);
  my $out=Bio::SeqIO->new(-format=>"genbank");
  my $i=0;
  while(my $seq=$in->next_seq){
    next if($seq->length < $$settings{min_length});
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
  --min_length  1   Minimum length of a contig
  "
}
