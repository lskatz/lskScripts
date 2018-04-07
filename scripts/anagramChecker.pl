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
  if($$settings{help} || @ARGV < 2){
    die "Usage: $0 referenceWord queryWord [queryWord2...]
    
    This script checks for anagrams as compared to the reference word.
    However, it will report false positives if there differences in
    letter counts.";
  }

  my $bitwiseLetters = bitwiseLetters($settings);

  # Simplify to uppercase for all
  @ARGV = map{uc($_)} @ARGV;

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

  return 0;
}

sub bitwiseLetters{
  my($settings)=@_;

  # Make the resolution as high as 10 letters per word.
  my $numLetters=26 * 10;

  my %bitwiseLetter;
  my $ordOffset=ord("A");
  for(my $i=0;$i<$numLetters;$i++){
    # Mod to find the letter of the alphabet
    my $mod = $i % 26;
    # String multiplier, e.g., A x 3 = AAA
    my $multiplier = int($i / 26)+1;
    # The chr that corresponds to the letter extended by $multiplier.
    my $key=chr($mod + $ordOffset) x $multiplier;

    # Power of 2 to take advantage of binary
    $bitwiseLetter{$key} = 2 ** $i;
  }

  return \%bitwiseLetter;
}

sub wordToBitwise{
  my($word,$bitwiseLetters,$settings)=@_;
  
  my $bitwise=0;
  my $sortedLetters = join("",sort{$a cmp $b} split(//,$word));

  while($sortedLetters=~/((.)\2*)/g){
    $bitwise = $bitwise | $$bitwiseLetters{$1};
  }
    
  return $bitwise;
}
