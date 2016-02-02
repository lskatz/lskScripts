#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Config::Simple;
use Getopt::Long;
use File::Spec::Functions qw/rel2abs/;

#my $baseDir="/scicomp/home/gzu2/projects/wgsStandards/accuracyVsCoverage/manyCoverages";
my $baseDir=".";

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help reps=i min_coverage=i max_coverage=i));
  $$settings{reps}||=1;
  $$settings{min_coverage}||=1;
  $$settings{max_coverage}||=$$settings{min_coverage};
  die usage() if($$settings{help});
  $$settings{config}=shift(@ARGV);
  die "ERROR: need config file" if(!$$settings{config});

  # Get configuration from cfg file
  my %config;
  Config::Simple->import_from($$settings{config},\%config);
  my $fixed=fixConfigValue(\%config);
  %config=%$fixed;


  for(my $i=$$settings{min_coverage};$i<=$$settings{max_coverage};$i++){
    for(my $rep=0; $rep<$$settings{reps}; $rep++){
      my $coverage=($i*1);

      # make a simulation directory
      my $simDir="cov$i/rep$rep";
      $simDir=rel2abs($simDir);
      system("mkdir -pv $simDir");
      # copy all necessary files except config
      system("cp -rv $baseDir/IlluminaErrorProfiles $baseDir/reference $baseDir/lyve-set.dnd $simDir/");
      
      # generate a custom config file
      my %newConfig=%config;
      # update some paths
      $newConfig{'default.output_dir'}="$simDir/out";
      $newConfig{'default.coverage'}=$coverage;
      $newConfig{'default.treefile_path'}="$simDir/$newConfig{'default.treefile_path'}";
      $newConfig{'default.base_genome_path'}="$simDir/".$newConfig{'default.base_genome_path'};
      for(qw(error_model1 error_model2)){
        $newConfig{"default.$_"}="$simDir/".$newConfig{"default.$_"};
      }
      
      open(CFG,">","$simDir/TTR.cfg") or die "ERROR: could not open $simDir/TTR.cfg for writing: $!";
      while(my($key,$value)=each(%newConfig)){
        $key=~s/^default\.//;
        if(ref($value) eq 'ARRAY'){
          $value=join(",",@$value);
        }
        print CFG "$key = $value\n";
      }
      close CFG;
    }
  }
}

sub fixConfigValue{
  my($value)=@_;
  if(ref($value) eq 'HASH'){
    while(my($k,$v)=each(%$value)){
      $$value{$k}=fixConfigValue($v);
    }
  } elsif(ref($value) eq 'ARRAY'){
    for(my $i=0;$i<@$value;$i++){
      $$value[$i]=fixConfigValue($$value[$i]);
    }
  } else {
    $value=~s/#.*$//;
    $value=~s/^\s+|\s+$//g;
  }
  return $value;
}

sub usage{
  "Create a bunch of different treetoreads projects
  Usage: $0 treetoreads.cfg
  --reps         1  Number of repetitions
  --min_coverage 1
  --max_coverage 1
  "
}
