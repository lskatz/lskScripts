#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long qw/GetOptions/;
use File::Basename qw/basename/;
use File::Copy qw/mv/;
use File::Temp qw/tempdir tempfile/;
use Data::Dumper;

local $0=basename $0;

sub logmsg{print STDERR "$0: @_\n"}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s taxid=i)) or die $!;
  die usage() if(!@ARGV);
  die "ERROR: need taxid" if(!defined $$settings{taxid});
  $$settings{tempdir}||=tempdir("$0.XXXXXX", TMPDIR=>1, CLEANUP=>1);

  my $suffix = '|kraken:taxid|'.$$settings{taxid};
  my $suffixRegex = qr/(.*?)(\|kraken:taxid\|\d+)*$/;
  
  logmsg "Taxid: $$settings{taxid}";
  for my $fasta(@ARGV){
    die "ERROR: cannot find $fasta" if(!-e $fasta);
    logmsg $fasta;

    my($tempfasFh, $tempfas)=tempfile("XXXXXX", SUFFIX=>".fasta", DIR=> $$settings{tempdir});
    my $in=Bio::SeqIO->new(-file=>$fasta);
    my $out=Bio::SeqIO->new(-fh=>$tempfasFh,-format=>"fasta");

    while(my $seq=$in->next_seq){
      my $id=$seq->id;
      $id=~s/$suffixRegex/$1$suffix/;
      $seq->id($id);
      $seq->description(" ");
      $out->write_seq($seq);
    }

    $in->close; 
    $out->close;
    close($tempfasFh);
    mv($tempfas,$fasta);
  }

  return 0;
}

sub usage{
  "$0: Perform in-place editing of fasta files for Kraken
  usage: $0 --taxid=taxid *.fasta

  --taxid     The taxon ID from NCBI (required)
  "
}

