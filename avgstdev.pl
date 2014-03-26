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
use lib "$FindBin::Bin";
use Statistics::Descriptive;

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

  my $stat=Statistics::Descriptive::Full->new();
  $stat->add_data(\@data);

  # overall metrics
  my $sum=sum(@data);
  my $min=min(@data);
  my $max=max(@data);
  print "Total: $sum\n";
  # normal dist stuff
  printf("Average: %f +/- %f\n",$stat->mean(),$stat->standard_deviation());
  # weighted distribution
  #print $stat->quantile(1)."\n";
  printf("Median: %0.2f [%0.2f,%0.2f] [%0.1f-%0.1f]\n",$stat->median(),$stat->quantile(1),$stat->quantile(3),$min,$max);

  if(@badData && $$settings{warnings}){
    warn "Warning: these data points were not used because they do not look like numbers: \n=>".join("\n=>",@badData)."\n";
  }

  return 0;
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
