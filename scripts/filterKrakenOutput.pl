#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw/GetOptions/;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use Data::Dumper;

local $0=basename $0;

sub logmsg{print STDERR "$0: @_\n"}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s taxid=s)) or die $!;
  die usage() if(!@ARGV);
  die "ERROR: need taxid" if(!defined $$settings{taxid});
  $$settings{tempdir}||=tempdir("$0.XXXXXX", TMPDIR=>1, CLEANUP=>1);

  my @taxid = split(/,/, $$settings{taxid});
  my @fastq = @ARGV;

  my $regex = "^".join('$|^', @taxid)."\$";
  $regex = qr/$regex/;

  my %readid = ();
  while(<STDIN>){
    my(undef, $readid, $taxid) = split(/\t/, $_);
    if($taxid =~ $regex){
      $readid{$readid} = 1;
    }
  }
  
  for my $f(@fastq){
    logmsg "Reading $f and filtering for $$settings{taxid}";
    open(my $fh, "zcat $f | ") or die "ERROR: could not gunzip $f: $!";
    while(my $id = <$fh>){
      my $entry = $id;
      $entry.=<$fh> for(1..3);

      $id =~ s/^@|\s+.*$//g; # remove @ and anything after whitespace

      if($readid{$id}){
        print $entry;
      }
    }
    close $fh;
  }

  return 0;
}
  
  #    cat out.kraken | perl -MData::Dumper -lane 'BEGIN{open($fh, "zcat SE-le_S12_L001_R1_001.fastq.gz | ") or die $!; while(my $id=<$fh>){$entry=$id; for(1..3){$entry.=<$fh>;} $id=~s/\s.*//; $id=~s/^\@//; chomp($id); $entry{$id}=$entry; } } my(undef,$readid,$taxid)=@F; next if($taxid !~ /^561$|^562$|^83334$/); print $entry{$readid};' | grep . > R1.subset.fastq

sub usage{
  "$0: Filter for reads matching a given taxon using kraken raw results
  usage: $0 --taxid=taxid in.fastq.gz < kraken.out > out.fastq

  --taxid     The taxon ID from NCBI (required)
              Can be comma-separated
  "
}

