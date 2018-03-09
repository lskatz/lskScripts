#!/usr/bin/env perl
# Lee Katz
# Moves a symbolic link or more to a target directory

use strict;
use warnings FATAL=>'all';
use File::Copy qw(mv);
use Getopt::Long qw(GetOptions);
use File::Spec;
use Cwd qw/realpath getcwd/;
use File::Basename;

my $target;
GetOptions('target-directory=s' => \$target);
die "$0 -t target_dir symlink1 symlink2 symlink3\n" unless $target && -d $target;

my $origDir=getcwd;
for (@ARGV) {
    unless (-l $_) {
        warn "$_ is not a symlink\n";
        next;
    }
    my $filename=fileparse $_;
    my $absPath=realpath($_);
    chdir $target;
      my $relPath=File::Spec->abs2rel($absPath);
      symlink $relPath, $filename;
    chdir $origDir;
    unlink $_;
}
