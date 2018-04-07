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
  GetOptions($settings,qw(help test operation=s));
  die usage() if($$settings{help});
  $$settings{operation}||="\$next - \$cur";
  #$$settings{operation}=quotemeta($$settings{operation});
  
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

    # do the math
    my $answer=eval($$settings{operation});
    if($@){
      die "ERROR: $$settings{operation} resulted in a failure: $@";
    }
    print "$answer\n";
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
  "Calculates the difference between rows or a custom arithmetic 
  Usage: sort -n numbers.txt | $0 > difference.txt
  -o 'custom arithmetic'
    Variables: \$cur is the first row in the iteration
               \$next is the second row
    Example:   \$next - \$cur
  "
}

