#!/usr/bin/env perl
#Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help test));
  die usage() if($$settings{help});
  
  if($$settings{test}){
    test($settings);
  } else {
    printDifferences($settings);
  }
  return 0;
}

sub printDifferences{
  my($settings)=@_;
  
  my $cur=<>; chomp($cur);
  while(my $next=<>){
    chomp($next);
    print(($next-$cur)."\n");
    $cur=$next;
  }
}

sub test{
  my($settings)=@_;
  my $cmd=" echo -e '1\n5\n7\n33\n33\n33\n37' | $0";
  logmsg "COMMAND:\n====\n$cmd\n====";
  system($cmd);
}

sub usage{
  "Calculates the difference between rows
  Usage: sort -n numbers.txt | $0 > difference.txt
  "
}

