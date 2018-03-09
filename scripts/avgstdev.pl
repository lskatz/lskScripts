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
use FindBin;
use lib "$FindBin::Bin/lib";
use Statistics::Descriptive;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help decimals=i warnings format=s ignore=i comment=s logarithmic lt=s gt=s));
  die usage() if($$settings{help});
  my $commentChar=$$settings{comment};
    $commentChar='#' if(!defined($commentChar));
  my $ignoreLines=$$settings{ignore} || 0;
  $$settings{decimals}//=2;
  $$settings{format}||="raw";

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

  # Warn about weird looking data
  if(@badData && $$settings{warnings}){
    warn "Warning: these data points were not used because they do not look like numbers: \n=>".join("\n=>",@badData)."\n";
  }

  printResults({
    data=>\@data,
  },$settings);
  
  return 0;
}

sub printResults{
  my($r,$settings)=@_;
  my $stat=Statistics::Descriptive::Full->new();
  $stat->add_data($$r{data});
  # overall metrics
  my $sum=sum(@{$$r{data}});
  my $min=min(@{$$r{data}});
  my $max=max(@{$$r{data}});
  my $d=$$settings{decimals};

  if($$settings{format} =~/^raw$/i){

    print "Total: $sum\n";
    # normal dist stuff
    printf("Average: %0.${d}f +/- %0.${d}f\n",$stat->mean(),$stat->standard_deviation());
    # weighted distribution
    #print $stat->quantile(1)."\n";
    printf("Median: %0.${d}f [%0.${d}f,%0.${d}f] [%0.${d}f-%0.${d}f]\n",$stat->median(),$stat->quantile(1),$stat->quantile(3),$min,$max);
    printf("MAD: %0.${d}f\n",medianAbsoluteDeviation($$r{data},$settings));

  }elsif($$settings{format} =~/^nonparametric$/i){
    printf("%0.${d}f ± %0.${d}f\n[%0.${d}f - %0.${d}f]\n",$stat->median(),medianAbsoluteDeviation($$r{data},$settings),$min,$max);
  }elsif($$settings{format} =~/^parametric$/i){
    printf("%0.${d}f ± %0.${d}f\n[%0.${d}f - %0.${d}f]\n",$stat->mean(),$stat->standard_deviation(),$min,$max);
  }
    
}

# Median absolute deviation (MAD):
#   1. All the absolute distances from the median
#   2. Take the median of those absolute distances
sub medianAbsoluteDeviation{
  my($data,$settings)=@_;
  my $stat=Statistics::Descriptive::Full->new();
  $stat->add_data(@$data);
  my $median=$stat->median();
  my @deviation;
  for(my $i=0;$i<@$data;$i++){
    push(@deviation,abs($$data[$i] - $median));
  }
  my $stat2=Statistics::Descriptive::Full->new();
  $stat2->add_data(@deviation);
  return $stat2->median;
}

sub usage{
  "Gets mean and stdev from stdin. And now median, [Q1,Q3], [min-max], and median absolute distance.
  Usage: cat 'listofdata.txt' | $0
  Another example: samtools depth bowtieAssembly.bam | cut -f 3| $0

  -w     for warnings
  -c     '#' to ignore lines with this leading comment character. Default: '#'
  -i     number to discard the first few lines
  -l     if the data are logarithmic and should be transformed before calculation
  -gt 0  Filter numbers not greater than this number
  -lt 0  Filter numbers not less than this number
  -d  2  Decimal digits to show
  -f     Output format. Possible options: raw(default), parametric, nonparametric
         Parametric displays avg/stdev
         Nonparametric displays median/MAD
  ";
}
