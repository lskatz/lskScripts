#!/usr/bin/env perl
# Author:  Lee Katz <lkatz@cdc.gov>
# Count the number of nucleotides in a sequence file

use strict;
use warnings;
use Bio::Perl;
use Getopt::Long;

exit(main());

sub main{
  
  my $settings={};
  GetOptions($settings,qw(help discoverNts)) or die $!;

  die usage() if($$settings{help});

  my @seq=@ARGV;
  die "ERROR: no sequence files given!\n".usage() if(!@seq);
  # Make an array of nts so that the ordering doesn't change
  my @nt=ntArray(\@seq,$settings);
  my @ntWithN=(@nt,"N");

  # print the header
  print $_."\t" for("file",@ntWithN);
  print "\n";

  # Count the ATCG for each parameter
  for my $seq(@seq){
    count($seq,\@nt,\@ntWithN);
  }
  return 0;
}

sub ntArray{
  my($seq,$settings)=@_;
  my @nt=qw(A T C G);
  if($$settings{discoverNts}){
    # make a long string of nucleotides
    my $concatSeq="";
    for my $seqfile(@$seq){
      my $seqin=Bio::SeqIO->new(-file=>$seqfile);
      while(my $seqObj=$seqin->next_seq){
        $concatSeq.=uc($seqObj->seq);
      }
    }
    # Figure out which nts are present
    my %nt;
    my $seqLength=length($concatSeq);
    for(my $i=0;$i<$seqLength;$i++){
      $nt{substr($concatSeq,$i,1)}=1;
    }
    @nt=keys(%nt);
  }
  return @nt;
}

sub count{
  my($arg,$ntArr,$ntWithN)=@_;
  # Read the sequence into an object
  my $in=Bio::SeqIO->new(-file=>$arg); 

  my %num;
  while(my $seq=$in->next_seq){
    $seq=$seq->seq; 
    # get the counts of nts
    for my $nt(@$ntArr){
      $num{$nt}=($seq=~s/$nt//gi); 
    } 
    $num{"N"}=length($seq);
  } 
  print $arg."\t";
  print $num{$_}."\t" for(@$ntWithN); # print values
  print "\n"; # newline to make it pretty
  $in->close;
}

sub usage{
  "Counts the number of each nucleotide in a sequence file
  Usage: $0 file.fasta [file2.fasta ...]
  --discoverNts  Look at the sequencing file first to figure out which Nts are present in the first place (slower)
  "
}
