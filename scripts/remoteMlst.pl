#!/usr/bin/env perl 

use warnings;
use strict;
use Cwd qw/abs_path/;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename dirname/;
use File::Copy qw/mv/;
use File::Path qw/rmtree/;
use File::Temp qw/tempdir/;
use File::Which qw/which/;

use version 0.77;
our $VERSION = '0.1.1';

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(tempdir=s outdir=s scheme=s numcpus=i help)) or die $!;
  $$settings{numcpus} ||= 1;
  $$settings{tempdir} ||= tempdir("remoteMlst.XXXXXX", TMPDIR=>1, CLEANUP=>1);
  usage() if($$settings{help} || !@ARGV);

  # Required options
  my $scheme = $$settings{scheme} or die "ERROR: need --scheme";
  my $outdir = $$settings{outdir} or die "ERROR: need --outdir";
  mkdir($outdir);
  mkdir("$outdir/raw");

  # Check for necessary executables
  for my $exe(qw(fasterq-dump mlst stringMLST.py)){
    $$settings{$exe} = which($exe);
    die "ERROR: could not find $exe in your PATH" if(!$$settings{$exe});
  }
  my $mlsthelp = `$$settings{mlst} --help`;
  my $mlstdb;
  if($mlsthelp =~ m|(/.+?mlst.fa)|){
    $mlstdb = $1;
  } else {
    die "ERROR: could not find the mlst database in the mlst help output";
  }

  my $formattedDb = formatDb($mlstdb, $scheme, $outdir, $settings);

  # Main loop
  for my $sra_run(@ARGV){
    my $outfile = "$outdir/raw/$sra_run.mlst.tsv";
    logmsg "Running MLST on $sra_run and saving to $outfile";
    if(-e $outfile){
      logmsg "Skipping $sra_run because $outfile already exists";
      next;
    }
    logmsg "Downloading $sra_run";
    my $dir = downloadSra($sra_run,$settings);
    runMlst($dir, $outfile, $formattedDb, $settings);
    rmtree($dir);
  }

  return 0;
}

sub formatDb{
  my($mlstdb, $scheme, $outdir, $settings) = @_;
  my $formattedDb = "$outdir/stringMLST.$scheme.db";
  if(-e "$formattedDb/_35.txt"){
    logmsg "Will not recreate database $formattedDb because it already exists";
    return $formattedDb;
  }

  mkdir($formattedDb);

  my $indir = abs_path(dirname($mlstdb) . "/../pubmlst/$scheme");
  if(!-e $indir){
    die "ERROR: could not find the pubmlst directory at $indir\n  Did you supply the correct scheme?";
  }

  # Find which loci we are looking for and put it into the stringMLST config
  my $config = "$$settings{tempdir}/config.txt";
  my @locusFile = glob("$indir/*.tfa");
  open(my $fh, ">", $config) or die "ERROR: could not write to $config: $!";
  print $fh "[loci]\n";
  for my $locusFile(@locusFile){
    my $locus = basename($locusFile, ".tfa");
    print $fh "$locus  $locusFile\n";
  }
  my $profile = "$indir/$scheme.txt";
  if(!-e $profile){
    die "ERROR: could not find profile file at $profile";
  }
  print $fh "[profile]\n";
  print $fh "profile  $profile\n";
  close $fh;
  # the config has been created into $config

  #system("cat $config");

  system($$settings{'stringMLST.py'}. " --buildDB --config $config -k 35 -P '$formattedDb/' | sed 's/^/[stringMLST.py] /' >&2");
  die "ERROR running stringMLST.py --buildDB" if $?;

  # => produces files _35.txt _profile.txt _weight.txt

  # Make this directory read-only to help make it thread safe
  for my $file(glob("$formattedDb/*")){
    # read-only
    chmod(0644, $file);
  }
  ## read-execute on the directory itself
  #chmod(0755, $formattedDb);

  return $formattedDb;
}

sub downloadSra{
  my($acc, $settings) = @_;
  my $tempdir = $$settings{tempdir} . "/$acc";
  mkdir $tempdir;
  system($$settings{'fasterq-dump'}." $acc --threads $$settings{numcpus} --outdir $tempdir --split-files --skip-technical | sed 's/^/[fasterq-dump] /' >&2");
  die "ERROR with fasterq-dump" if $?;
  system("gzip -1v $tempdir/*.fastq 1>&2");
  die "ERROR with gzip" if $?;

  return $tempdir;
}

sub runMlst{
  my($indir, $out, $db, $settings) = @_;
  my $tmpout = $$settings{tempdir}."/mlst.out";
  my @fastq  = glob("$indir/*.fastq.gz");
  my $sample = basename($fastq[0], qw(_1.fastq.gz));
  if(-e $out){
    logmsg "Skipping $out because it already exists";
    return $out;
  }

  # If we have only one fastq file, then I'll need to rethink this.
  if(@fastq == 1){
      system($$settings{'stringMLST.py'} . " --predict -1 $fastq[0] -s -P $db/ -k 35 -o $tmpout | sed 's/^/[stringMLST.py] /' >&2");
  }
  if(@fastq > 2){
    logmsg "ERROR: expected 2 fastq files but found ".scalar(@fastq);
    return "";
  }
  if(@fastq == 2){
    system($$settings{'stringMLST.py'} . " --predict -1 $fastq[0] -2 $fastq[1] --paired -P $db/ -k 35 -o $tmpout | sed 's/^/[stringMLST.py] /' >&2");
  }
  
  # => produces, e.g.,
  # Sample  dnaE    dtdS    gyrB    pntA    pyrC    recA    tnaA    ST
  # SRR19910447     3       4       4       29      4       19      22      3

  mv($tmpout, $out) or die "ERROR: could not move results to $out: $!";

  return $out;
}

sub usage{
  print "$0: runs MLST on an NCBI SRA run accession and returns the sequence type per sample
  Usage: $0 [options] acc1 acc2 acc3 ...
  --scheme  The MLST scheme to use.
            Run 'mlst --list' to see available schemes.
  --outdir  Output directory
  --numcpus Number of CPUs to use. Default: 1
  --help    This useful help menu
  \n";
  exit 0;
}
