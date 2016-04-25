#!/usr/bin/env perl

# validates a fastq file.  Is it in the right format?
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/basename/;

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
  my $readLength=0;
  while(<>){
    $i++;
    my $mod=$i%4;

    # ID line
    if($mod==1){
      if(substr($_,0,1) ne '@'){
        return errorMsg($i,$_,$settings);
      }
    } 
    # seq line
    elsif($mod==2){
      chomp;
      if($_=~/[^A-Z\.]/i){
        return errorMsg($i,$_,$settings);
      }
      $readLength=length($_);
    } 
    # plus line
    elsif($mod==3){
      if(substr($_,0,1) ne '+'){
        return errorMsg($i,$_,$settings);
      }
    } 
    # qual line
    elsif($mod==0){
      chomp or return errorMsg($i,$_,$settings);
      # This usually captures incomplete downloads.
      if($readLength != length($_)){
        return "Seq length is not the same as qual length ".errorMsg($i,$_,$settings);
      }
    }
  }

  # Did we get four lines per read?
  if($i % 4 > 0){
    return "The last read is incomplete ". errorMsg($i,"The last line",$settings);
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
  $0=basename $0;
  "Validates the format of a fastq file and returns an appropriate exit code.
  The fastq file must be in a 4-line-per-read format.

  Usage: $0 < reads.fastq
  -h for help
  -v verbose
  Examples:
    gunzip -c reads.fastq.gz | $0
    $0 < reads.fastq
    $0 -v < reads.fastq # to print where it went wrong
    "
}
