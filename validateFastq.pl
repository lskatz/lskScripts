#!/usr/bin/env perl

# validates a fastq file.  Is it in the right format?
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long qw/GetOptions/;
use File::Basename qw/basename/;
use List::Util qw/sum/;
require v5.12;

sub logmsg{local $0=basename($0); print STDERR "$0: @_\n";}

exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help verbose pe|paired-end|pairedend min-quality|min-qual=i min-length=i));
  $$settings{pe}//=0;
  $$settings{'min-quality'}//=0;
  $$settings{'min-length'}//=0;
  die usage() if($$settings{help});

  my $err=validate($settings);
  if($err){
    logmsg $err;
    return 1;
  }elsif($$settings{verbose}){
    logmsg "This fastq file seems to be intact.";
  }

  return 0;
}

# Returns "" if everything's good.
# Returns a string if it's not.
sub validate{
  my($settings)=@_;

  my $linesPerEntry=4;
  if($$settings{pe}){
    $linesPerEntry=8;
  }

  # Check if any of the detailed/slow options are checked
  # so that only one if statement has to be used in the loop
  # collectively.
  my $fastqEntryOptions=0;
  if($$settings{'min-length'} || $$settings{'min-quality'}){
    $fastqEntryOptions=1;
  }
  my $minLength=$$settings{'min-length'};
  my $minQual=$$settings{'min-quality'};

  # validate the fastq file
  my $i=0;
  my $readLength=0;
  my $qualLine="";
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
      $qualLine=$_;

      if($fastqEntryOptions==1){
        my $err=checkDetailedOptions($readLength,$qualLine,$minQual,$minLength,$settings);
        return $err.errorMsg($i,$_,$settings) if($err);
      }
    }


  }

  # Did we get four lines per read?
  if($i % 4 > 0){
    return "The last fastq entry does not have four lines.";
  }

  if($i % $linesPerEntry > 0){
    return "The last fastq entry is unpaired.";
  }
  if($i == 0){
    return "There are zero reads!";
  }

  return "";
}

# Any slow/detailed fastq entry error check goes here.
sub checkDetailedOptions{
  my($readLength,$qualLine,$minQual,$minLength,$settings)=@_;
  if($readLength < $minLength){
    return "Seq length is less than the min length\n";
  }
  
  if($minQual > 0){
    my @qual;
    for my $q (split(//,$qualLine)){
      push(@qual,ord($q));
    }
    my $avgQual=sum(@qual)/$readLength - 33;
    if($avgQual < $minQual){
      $avgQual=sprintf("%0.1f",$avgQual);
      return "A read's quality ($avgQual) is less than the min-quality ($minQual).\n";
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
  $0=basename $0;
  "Validates the format of a fastq file and returns an appropriate exit code.
  The fastq file must be in a 4-line-per-read format.

  Usage: $0 < reads.fastq
  --help        for help
  --verbose     verbose
  --pe          Check for interleaved paired end

  Options that might slow it down
  --min-length  1
  --min-quality 1

  Examples:
    zcat reads.fastq.gz | $0
    $0 < reads.fastq
    "
}
