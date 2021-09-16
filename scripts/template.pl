#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  usage() if($$settings{help});

  return 0;
}

sub usage{
  print "$0: does something
  Usage: $0 [options] arg1
  --help   This useful help menu
  ";
  exit 0;
}
