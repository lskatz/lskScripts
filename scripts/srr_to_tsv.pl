#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use File::Copy qw/mv/;
use File::Which qw/which/;

use threads;
use Thread::Queue;

use version 0.77;
our $VERSION = '0.1.1';

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(numcpus=i tmpdir=s help)) or die $!;
  usage() if($$settings{help});

  if(-t STDIN){
    logmsg "ERROR: no STDIN detected";
    usage();
  }

  if(!which("esearch")){
    die "ERROR: could not find esearch in PATH";
  }

  $$settings{numcpus} ||= 1;
  $$settings{tmpdir}  ||= tempdir("$0.XXXXXX", TMPDIR=>1, CLEANUP=>1);
  mkdir($$settings{tmpdir});

  convertSrrToTsv($settings);

  return 0;
}

sub convertSrrToTsv{
  my($settings) = @_;

  my @SRR;
  while(my $srrLine = <STDIN>){
    chomp($srrLine);
    for my $srr(split(/\s+/, $srrLine)){
      push(@SRR, $srr);
    }
  };
  close STDIN;

  my $Q = Thread::Queue->new(@SRR);
  my $printerQ = Thread::Queue->new;
  my $printer = threads->new(\&printer, $printerQ, $settings);

  my @thr;
  for(my $i=0;$i<$$settings{numcpus};$i++){
    $thr[$i] = threads->create(\&convertWorker, $Q, $printerQ, $settings);
  }

  # Terminate threads
  for(@thr){
    $Q->enqueue(undef);
  }
  for(@thr){
    $_->join;
  }
  $printerQ->enqueue(undef);
  $printer->join;

  return scalar(@SRR);
}

sub printer{
  my($Q, $settings) = @_;
  while(defined(my $line = $Q->dequeue)){
    print $line . "\n";
  }
}

sub convertWorker{
  my($Q, $printer, $settings) = @_;

  while(defined(my $srr = $Q->dequeue)){
    my $info = convertSrr($srr, $settings);
    my $line = join("\t",
      $$info{srr},
      $$info{biosample},
      $$info{strain},
    );

    $printer->enqueue($line);
  }
  
  return 1;
}

sub convertSrr{
  my($srr, $settings) = @_;

  my $tempdir = "$$settings{tmpdir}/$srr";
  mkdir($tempdir);

  my %info;
  $info{srr} = $srr;
  $info{tempdir} = $tempdir;
  
  # Create the edirect XML files into the temp dir.
  # First into .xml.tmp files and then into .xml files
  # to help with checking whether the files exist and are intact.
  my @esearchStat = stat("$tempdir/esearch.xml.gz");
  my @biosampleStat = stat("$tempdir/biosample.xml.gz");
  #print Dumper \@biosampleStat, $biosampleStat[7];
  if(!-e "$tempdir/esearch.xml.gz" || $esearchStat[7] < 100){
    system("esearch -db sra -query '$srr' | gzip -c > $tempdir/esearch.xml.gz.tmp");
    die if $?;
    # TODO show a warning if "count" is zero in the xml
    mv("$tempdir/esearch.xml.gz.tmp", "$tempdir/esearch.xml.gz");
  }
  if(!-e "$tempdir/biosample.xml.gz" || $biosampleStat[7] < 100){
    system("zcat $tempdir/esearch.xml.gz | elink -target biosample | efetch -format xml | xtract -format | gzip -c > $tempdir/biosample.xml.gz.tmp");
    die if $?;
    mv("$tempdir/biosample.xml.gz.tmp", "$tempdir/biosample.xml.gz");
  }

  $info{biosample}   = qx(zcat $tempdir/biosample.xml.gz | xtract -pattern BioSample -element BioSample\@accession);
  chomp($info{biosample});
  if(!$info{biosample}){
    $info{biosample} = "MISSING";
  }

  $info{strain} = qx(zcat $tempdir/biosample.xml.gz | grep -Ev 'is_primary|SRA' | xtract -pattern BioSample -block Ids -if Id\@db_label -equals "Sample name" -element Id);
  # If we don't have the strain name then it might be in a different label like Id@db=CFSAN.
  # If so, check for whatever ID is there.
  if(!$info{strain}){
    $info{strain} = qx(zcat $tempdir/biosample.xml.gz | grep -Ev 'is_primary|SRA' | xtract -pattern BioSample -block Ids -element Id);
    # If there is still no ID, label as MISSING.
    if(!$info{strain}){
      $info{strain} = "MISSING";
    }
  }
  # If there are somehow multiple IDs, just get the first.
  $info{strain} = (split(/\s+/, $info{strain}))[0];
  chomp($info{strain});

  return \%info;

}

sub usage{
  print "$0: converts a set of SRRs into a spreadsheet of accessions
  Accessions are separated by any whitespace including spaces or newlines.
  Usage: $0 [options] < SRRs.txt
         cat SRRs.txt | $0 [options]
  --tmpdir  A temporary directory
  --help     This useful help menu
  \n";
  exit 0;
}
