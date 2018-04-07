#!/usr/bin/env perl

# Give the approximate lenght of time it took to make all files in a directory

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/basename/;
use File::Find qw/find/;

local $0=basename $0;

my $settings={};
GetOptions($settings,qw(help exclude=s include=s verbose)) or die $!;
die usage() if(!@ARGV || $$settings{help});
my $exclude=$$settings{exclude}||0;
my $include=$$settings{include}||0;

for my $dir(@ARGV){
  die "ERROR: $dir is not a directory" if(!-d $dir);

  my $oldest= ~0;
  my $newest=0;
  find({no_chdir=>1, wanted=>sub{
    return if(-d $File::Find::name);
    if($exclude && $File::Find::name =~ /$exclude/){
      print "Excluding: $File::Find::name\n" if($$settings{verbose});
      return;
    }
    if($include && $File::Find::name !~ /$include/){
      print "Not including: $File::Find::name\n" if($$settings{verbose});
      return;
    }
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
  Gives the number of seconds between the oldest and 
  newest files

  --verbose
  --exclude  PATTERN  Supply a regex pattern to ignore
                      certain filenames.
  "
}
