#!/usr/bin/env perl

use warnings;
use strict;
use File::Path;
use Bio::Perl;
use File::Basename;

my $out=Bio::SeqIO->new(-format=>"fasta");
my @file=@ARGV or die("need files to convert\n");
foreach my $f (@file){
  my $in=Bio::SeqIO->new(-format=>"lasergene",-file=>$f);
  while(my $seq=$in->next_seq){
    my($id,$desc);
    ($id)=fileparse($f);
    ($id,$desc)=split(/\s+/,$id,2);

    $seq->id($id);
    $seq->desc($desc);
    $out->write_seq($seq);
  }
}

sub trim{
  my $str=shift;
  $str=~s/^\s+|\s+$//g;
  return $str;
}
