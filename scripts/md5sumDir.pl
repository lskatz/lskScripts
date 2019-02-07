#!/usr/bin/env perl
use strict;
use warnings;
use Digest::MD5 qw/md5_hex/;
use File::Slurp qw/read_file/;
use threads;

my $numcpus = 24;

my @file;
for my $dir(@ARGV){
  push(@file, 
    `find $dir -type f`
  );
}
chomp(@file);

my $num_files_per_thread = int(scalar(@file)/$numcpus) + 1;
my @thr;
for my $i(0..$numcpus - 1){
  my @subfile = splice(@file, 0, $num_files_per_thread);
  $thr[$i] = threads->new(sub{
    my @hex;
    for my $file(@subfile){
      my $content = read_file($file);
      push(@hex, md5_hex($content));
    }
    return \@hex;
  });
}

my @hex;
for my $thr(@thr){
  my $hexSubArr = $thr->join;
  push(@hex, @$hexSubArr);
}
my $finalHex = md5_hex(join("\n", sort @hex)."\n");
print $finalHex."\n";
