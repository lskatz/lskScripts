#!/usr/bin/env perl

# validates a fastq file.  Is it in the right format?
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;

exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help verbose));
  $$settings{lineNum}=0;
  die usage() if($$settings{help});

  my $err=validate($settings);
  if($err){
    print "ERROR! ".$err."\n";
    return 1;
  }
  return 0;
}

# Returns "" if everything's good.
# Returns a string if it's not.
sub validate{
  my($settings)=@_;
  # validate the fastq file
  my $i=0;
  while(<>){
    $i++;
    my $mod=$i%4;
    if($mod==1){
      if(substr($_,0,1) ne '@'){
        return errorMsg($i,$_,$settings);
      }
    } elsif($mod==2){
      chomp;
      if($_=~/[^A-Z\.]/i){
        return errorMsg($i,$_,$settings);
      }
    } elsif($mod==3){
      if(substr($_,0,1) ne '+'){
        return errorMsg($i,$_,$settings);
      }
    }
  }
  return "";
}

sub errorMsg{
  my($lineNum,$line,$settings)=@_;
  my $error="Input line number: $lineNum
  offending line is:
  $line";
  return $error;
}

sub usage{
  $0=`basename $0`; chomp($0);
  "Validates the format of a fastq file and returns an appropriate exit code. The fastq file must be in a 4-line-per-read format.
  Usage: $0 < reads.fastq
  -h for help
  -v verbose
  Examples:
    gunzip -c reads.fastq.gz | $0
    $0 < reads.fastq
    $0 -v < reads.fastq # to print where it went wrong
    "
}
