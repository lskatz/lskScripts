#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my %length;
my %percentage;
open(my $fh, $ARGV[0]) or die "ERROR: could not read $ARGV[0]: $!";
while(<$fh>){
  chomp;
  my($classified,$seqname,$taxid,$length,$kmerTaxid)=split(/\t/,$_);
  if($classified eq 'U'){
    $percentage{'unclassified'}+=$length;
  } else {
    $length{$seqname}=$length;
  }
}
close $fh;


# kraken-translate but tally all the sequence lengths
open(my $translateFh, "kraken-translate $ARGV[0] | ") or die "ERROR: could not run kraken-translate on $ARGV[0]:$!";
while(<$translateFh>){
  chomp;
  my($seqname,$taxonomyString)=split(/\t/,$_);
  $taxonomyString=~s/\s+/_/g;
  $taxonomyString=~s/;/\t/g;
  $percentage{$taxonomyString}+=$length{$seqname};
}
close $translateFh;

# Make the file
while(my($taxonomyString,$sliceOfPie)=each(%percentage)){
  print join("\t",$sliceOfPie,$taxonomyString)."\n";
}
