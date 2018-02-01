#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw/GetOptions/;
use File::Basename qw/basename/;
use Data::Dumper;

local $0 = basename $0;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  if($$settings{help}){
    die "Usage: $0 referenceWord queryWord [queryWord2...]
    
    This script checks for anagrams as compared to the reference word.
    However, it will report false positives if there differences in
    letter counts.";
  }

  my $bitwiseLetters = bitwiseLetters($settings);

  my $refWord=shift(@ARGV);
  my $refBitwise = wordToBitwise($refWord,$bitwiseLetters,$settings);
  for my $queryWord(@ARGV){
    my $queryBitwise = wordToBitwise($queryWord,$bitwiseLetters,$settings);
    if($refBitwise == $queryBitwise){
      print "$queryWord is an anagram of $refWord\n";
    } else {
      print "$queryWord is not an anagram of $refWord\n";
    }
  }
}

sub bitwiseLetters{
  my($settings)=@_;

  my %bitwiseLetter;
  my $ordOffsetUc=ord("A");
  my $ordOffsetLc=ord("a");
  for(my $i=0;$i<26;$i++){
    $bitwiseLetter{chr($i + $ordOffsetUc)} = 2 ** $i;
    $bitwiseLetter{chr($i + $ordOffsetLc)} = 2 ** $i;
  }

  return \%bitwiseLetter;
}

sub wordToBitwise{
  my($word,$bitwiseLetters,$settings)=@_;
  
  my $bitwise=0;
  for my $letter(split(//,$word)){
    $bitwise = $bitwise | $$bitwiseLetters{$letter};
  }
  return $bitwise;
}
