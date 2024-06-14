#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use File::Copy qw/cp/;
use File::Temp qw/tempdir/;

use version 0.77;
our $VERSION = '0.1.1';

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help db=s numcpus=i)) or die $!;
  usage() if($$settings{help});

  $$settings{tempdir} //= tempdir("kaptive.XXXXXX", TMPDIR=>1, CLEANUP=>1);
  $$settings{db} or die "ERROR: need --db set to the database with the kaptive database";
  $$settings{numcpus} ||= 1;

  for my $fasta(@ARGV){
    my $outdir = runKaptive($fasta, $$settings{db}, $settings);
    
    # Get the results
    open(my $kFh, "<", "$outdir/k.log") or die "ERROR: could not read $outdir/k.log: $!";
    my $k = <$kFh>;
    chomp($k);
    close $kFh;
    open(my $oFh, "<", "$outdir/o.log") or die "ERROR: could not read $outdir/o.log: $!";
    my $o = <$oFh>;
    chomp($o);
    close $oFh;

    # Parse the results
    $k =~ s/.*?K/K/;
    $o =~ s/.*?O/O/;
    print join("\t", basename($fasta), $k, $o)."\n";
  }

  return 0;
}

sub runKaptive{
  my($fasta, $db, $settings) = @_;
  logmsg "Running Kaptive on $fasta";
  
  # Staging area
  my $tmpdir = "$$settings{tempdir}/".basename($fasta, qw(.fasta .fa));
  mkdir $tmpdir;

  # Input and output directories
  my $tmpdirIn  = "$tmpdir/in";
  my $tmpdirOut = "$tmpdir/out";
  mkdir $tmpdirIn;
  mkdir $tmpdirOut;

  my $tmpfasta = "$tmpdirIn/in.fasta";
  cp($fasta, $tmpfasta) or die "ERROR: could not copy $fasta to $tmpfasta: $!";

  my $oGbk = "$db/VibrioPara_Kaptivedb_O.gbk";
  my $kGbk = "$db/VibrioPara_Kaptivedb_K.gbk";

  system("kaptive.py --threads $$settings{numcpus} -k $oGbk -a $tmpfasta -o $tmpdirOut/o > $tmpdirOut/o.log");
  die if $?;
  system("kaptive.py --threads $$settings{numcpus} -k $kGbk -a $tmpfasta -o $tmpdirOut/k > $tmpdirOut/k.log");
  die if $?;

  return $tmpdirOut;
}

sub usage{
  print "$0: runs Kaptive on a set of fasta files
  Usage: $0 [options] *.fasta > out.tsv
  --db      Database directory for Kaptive containing *.gbk
  --numcpus Number of threads to use (default: 1)
  --help    This useful help menu
  \n";
  exit 0;
}
