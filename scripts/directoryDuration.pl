#!/usr/bin/env perl

# Give the approximate lenght of time it took to make all files in a directory

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/basename/;
use File::Find qw/find/;

local $0=basename $0;

my $settings={};
GetOptions($settings,qw(help verbose)) or die $!;
die usage() if(!@ARGV || $$settings{help});


for my $dir(@ARGV){
  die "ERROR: $dir is not a directory" if(!-d $dir);

  my $oldest= ~0;
  my $newest=0;
  find({no_chdir=>1, wanted=>sub{
    return if(-d $File::Find::name);
    my ($dev,$ino,$mode,$nlink,$uid,$gid,$rdev,$size,
            $atime,$mtime,$ctime,$blksize,$blocks)
      =stat($File::Find::name);

   if(!$mtime){
     print "WARNING: no timestamp $File::Find::name\n" if($$settings{verbose});
     return;
   }
   if($mtime < $oldest){
     print "oldest $File::Find::name\n" if($$settings{verbose});
     $oldest=$mtime;
   }
   if($mtime > $newest){
     print "newest $File::Find::name\n" if($$settings{verbose});
     $newest=$mtime;
   }
  }},$dir);

  my $duration=$newest-$oldest;
  print "$dir\t$duration\n";
}


exit 0;

sub usage{
  "Usage: $0 dir
  --verbose
  Gives the number of seconds between the oldest and newest files
  "
}
