#!/usr/bin/env perl

# Author: Lee Katz <lkatz@cdc.gov>
# calculate the standard deviation and mean from stdin
# TODO calculate median absolute distance as well, for the median

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw/sum min max/;
use Scalar::Util qw/looks_like_number/;
#use Statistics::Descriptive;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help warnings ignore=i comment=s logarithmic lt=s gt=s));
  die usage() if($$settings{help});
  my $commentChar=$$settings{comment};
    $commentChar='#' if(!defined($commentChar));
  my $ignoreLines=$$settings{ignore} || 0;

  my @data;
  my @badData;
  my $i=0;
  while(<>){
    # skip $i number of lines
    next if($i++<$ignoreLines);
    s/^\s*|\s*$//g; # trim whitespace
    # skip lines with comments (after considering $i number of lines)
    next if($commentChar && /^$commentChar/);
    # each possible number is whitespace delimited, so split on whitespace
    my @tmp=split(/\s+/);
    for(@tmp){
      s/\s+//g;
      next if(/^$/);

      # categorize the data and remove non-numbers
      if(looks_like_number($_)){
        next if(defined($$settings{gt}) && $_<=$$settings{gt});
        next if(defined($$settings{lt}) && $_>=$$settings{lt});
        $_=exp($_) if($$settings{logarithmic});
        push(@data,$_);
      } else {
        push(@badData,$_);
      }
    }
  }

  my($mean,$std)=meanStdev(\@data,$settings);
  my $sum=sum(@data);

  print "$mean +/- $std\n";
  print "total: $sum\n";

  my $q1=median(\@data,0.25);
  my $q3=median(\@data,0.75);
  print "median: ". median(\@data)." [$q1,$q3]\n";
  print "Range: [".min(@data).'-'.max(@data)."]\n";

  if(@badData && $$settings{warnings}){
    warn "Warning: these data points were not used because they do not look like numbers: \n=>".join("\n=>",@badData)."\n";
  }

  return 0;
}

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}


sub meanStdev{
  my($data,$settings)=@_;
  my $avg=average($data);
  my $stdev=stdev($data);
  return($avg,$stdev);
}

# https://github.com/fhcrc/bdes/blob/master/Perl/Median.pl
sub median(){
    my($array,$quartile)=@_;
    $quartile||=0.5;
    my @array = sort {$a <=> $b } @$array;
    my $size = scalar(@array);
    my $i=$quartile * ($size+1);

    # if a decimal, then average above and below
    my $value;
    if($i!=int($i)){
      $value=(($array[$i]+$array[$i-1])/2)
    } else {
      $value=$array[$i];
    }
    return $value;
}

sub usage{
  "Gets mean and stdev from stdin. And now median.
  Usage: cat 'listofdata.txt' | $0
  Another example: samtools depth bowtieAssembly.bam | cut -f 3| $0

  -w for warnings
  -c '#' to ignore lines with this leading comment character. Default: '#'
  -i number to discard the first few lines
  -l if the data are logarithmic and should be transformed before calculation
  -gt 0 Filter numbers not greater than this number
  -lt 0 Filter numbers not less than this number
  ";
}
