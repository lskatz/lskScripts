#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/dirname basename/;
use File::Temp qw/tempdir tempfile/;
use Getopt::Long;
use File::Spec;

use Bio::SeqIO;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s)) or die $!;
  $$settings{tempdir}||=tempdir("$0.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  mkdir($$settings{tempdir}) if(!-d $$settings{tempdir});

  my($intable,@fasta)=@ARGV;
  my $numFasta=@fasta;

  die usage() if(!@fasta);
  
  # Find Primer3 and also the Tm parameters
  $$settings{primer3}=`which primer3_core 2>/dev/null` || `which primer3 2>/dev/null` || die "ERROR: could not find primer3_core or primer3 in path";
  chomp($$settings{primer3});
  my $exeDir=dirname($$settings{primer3});
  my $tmParamsPath="$exeDir/../test/primer3_config";
  if(!-e $tmParamsPath){
    logmsg "Could not find $tmParamsPath";
    $tmParamsPath="$exeDir/../src/primer3_config";
  }
  if(!-e $tmParamsPath){
    logmsg "Could not find $tmParamsPath";
    die "ERROR: could not find primer3_config! Should be in in the executable folder under src or test.";
  }
  $tmParamsPath=File::Spec->rel2abs($tmParamsPath);
  $$settings{tmParamsPath}="$tmParamsPath/";
  
  # Read in the SNP positions
  open(my $fh, $intable) or die "ERROR: could not read $intable: $!";
  # Get genome names
  my $header=<$fh>;
  chomp($header);
  my (undef, undef, undef, @genome)=split(/\t/, $header);
  my $numGenomes=@genome;

  # Loop through each genomic site
  while(<$fh>){
    chomp;
    my($seqname, $pos)=split(/\t/, $_);

    my @siteConfig;
    for(my $i=0;$i<$numFasta;$i++){
      my $confFile=primer3Config($fasta[$i], $seqname, $pos, $settings);
      push(@siteConfig,$confFile);
    }

    # Run Primer3 on all genomes at the same site
    primer3(\@siteConfig, $settings);

  }

  return 0;
}

sub primer3{
  my($configArray, $settings)=@_;

  my $primer3Dir=tempdir("primer3Out.XXXXXX", DIR=>$$settings{tempdir});
  my $allConfig="$primer3Dir/primer3.conf";

  # Concat all config files
  open(my $allConfigFh, ">", $allConfig) or die "ERROR: could not write to $allConfig: $!";
  for my $c(@$configArray){
    open(my $cFh, "<", $c) or die "ERROR: could not read $c: $!";
    while(<$cFh>){
      print $allConfigFh $_;
    }
    close $cFh;
  }
  close $allConfigFh;
  die;

  system("cd $primer3Dir && $$settings{primer3} < $allConfig > $primer3Dir/log");
  die "ERROR: with $$settings{primer3}" if $?;

}


sub primer3Config{
  my($fasta, $seqname, $pos, $settings)=@_;

  # Make a file with the parameters
  my($paramFh, $paramFile)=tempfile("primer3Config.XXXXXX",DIR=>$$settings{tempdir}, SUFFIX=>".conf");
  my $param="$$settings{tempdir}/".basename($fasta).".txt";

  # Get a subseq of the fasta. 500bp?
  my $start=$pos-499;
  my $end  =$pos+500;

  my $inseq=Bio::SeqIO->new(-file=>$fasta);
  my $subseq="";
  while(my $seq=$inseq->next_seq){
    if($seq->id eq $seqname){
      $subseq=$seq->subseq($start, $end);
      last;
    }
  }
  $inseq->close;

  # Make the params file
  my %param=(
    SEQUENCE_ID                       => join(":", basename($fasta,qw(.fasta)),$seqname,$start,$end),
    SEQUENCE_TEMPLATE                 => $subseq,
    SEQUENCE_TARGET                   => "500,1",
    PRIMER_TASK                       => "pick_sequencing_primers",
    PRIMER_PICK_LEFT_PRIMER           => 1,
    PRIMER_PICK_INTERNAL_OLIGO        => 0,
    PRIMER_PICK_RIGHT_PRIMER          => 1,
    PRIMER_OPT_SIZE                   => 20,
    PRIMER_MIN_SIZE                   => 18,
    PRIMER_MAX_SIZE                   => 23,
    PRIMER_MAX_NS_ACCEPTED            => 1, 
    PRIMER_PRODUCT_SIZE_RANGE         => "75-500",
    P3_FILE_FLAG                      => 1,
    SEQUENCE_INTERNAL_EXCLUDED_REGION => "400,200",
    PRIMER_EXPLAIN_FLAG               => 1,
    PRIMER_THERMODYNAMIC_PARAMETERS_PATH => $$settings{tmParamsPath},
  );

  # web params that I obtained on March 14, 2017
  my @webParam=qw(
  PRIMER_FIRST_BASE_INDEX=1
  PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1
  PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0
  PRIMER_PICK_LEFT_PRIMER=1
  PRIMER_PICK_INTERNAL_OLIGO=0
  PRIMER_PICK_RIGHT_PRIMER=1
  PRIMER_LIBERAL_BASE=1
  PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0
  PRIMER_LOWERCASE_MASKING=0
  PRIMER_PICK_ANYWAY=1
  PRIMER_EXPLAIN_FLAG=1
  PRIMER_TASK=pick_sequencing_primers
  PRIMER_MIN_QUALITY=0
  PRIMER_MIN_END_QUALITY=0
  PRIMER_QUALITY_RANGE_MIN=0
  PRIMER_QUALITY_RANGE_MAX=100
  PRIMER_MIN_SIZE=18
  PRIMER_OPT_SIZE=20
  PRIMER_MAX_SIZE=23
  PRIMER_MIN_TM=57.0
  PRIMER_OPT_TM=59.0
  PRIMER_MAX_TM=62.0
  PRIMER_PAIR_MAX_DIFF_TM=5.0
  PRIMER_TM_FORMULA=1
  PRIMER_PRODUCT_MIN_TM=-1000000.0
  PRIMER_PRODUCT_OPT_TM=0.0
  PRIMER_PRODUCT_MAX_TM=1000000.0
  PRIMER_MIN_GC=30.0
  PRIMER_OPT_GC_PERCENT=50.0
  PRIMER_MAX_GC=70.0
  PRIMER_PRODUCT_SIZE_RANGE=75-500
  PRIMER_NUM_RETURN=10
  PRIMER_MAX_END_STABILITY=9.0
  PRIMER_MAX_LIBRARY_MISPRIMING=12.00
  PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=20.00
  PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00
  PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00
  PRIMER_MAX_SELF_ANY_TH=45.0
  PRIMER_MAX_SELF_END_TH=35.0
  PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0
  PRIMER_PAIR_MAX_COMPL_END_TH=35.0
  PRIMER_MAX_HAIRPIN_TH=24.0
  PRIMER_MAX_TEMPLATE_MISPRIMING=-1.0
  PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=-1.0
  PRIMER_MAX_SELF_ANY=8.00
  PRIMER_MAX_SELF_END=3.00
  PRIMER_PAIR_MAX_COMPL_ANY=8.00
  PRIMER_PAIR_MAX_COMPL_END=3.00
  PRIMER_MAX_NS_ACCEPTED=0
  PRIMER_MAX_POLY_X=4
  PRIMER_INSIDE_PENALTY=-1.0
  PRIMER_OUTSIDE_PENALTY=0
  PRIMER_GC_CLAMP=0
  PRIMER_MAX_END_GC=5
  PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=3
  PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3
  PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION=7
  PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4
  PRIMER_SALT_MONOVALENT=50.0
  PRIMER_SALT_CORRECTIONS=1
  PRIMER_SALT_DIVALENT=1.5
  PRIMER_DNTP_CONC=0.6
  PRIMER_DNA_CONC=50.0
  PRIMER_SEQUENCING_SPACING=500
  PRIMER_SEQUENCING_INTERVAL=250
  PRIMER_SEQUENCING_LEAD=50
  PRIMER_SEQUENCING_ACCURACY=20
  PRIMER_WT_SIZE_LT=1.0
  PRIMER_WT_SIZE_GT=1.0
  PRIMER_WT_TM_LT=1.0
  PRIMER_WT_TM_GT=1.0
  PRIMER_WT_GC_PERCENT_LT=0.0
  PRIMER_WT_GC_PERCENT_GT=0.0
  PRIMER_WT_SELF_ANY_TH=0.0
  PRIMER_WT_SELF_END_TH=0.0
  PRIMER_WT_HAIRPIN_TH=0.0
  PRIMER_WT_TEMPLATE_MISPRIMING_TH=0.0
  PRIMER_WT_SELF_ANY=0.0
  PRIMER_WT_SELF_END=0.0
  PRIMER_WT_TEMPLATE_MISPRIMING=0.0
  PRIMER_WT_NUM_NS=0.0
  PRIMER_WT_LIBRARY_MISPRIMING=0.0
  PRIMER_WT_SEQ_QUAL=0.0
  PRIMER_WT_END_QUAL=0.0
  PRIMER_WT_POS_PENALTY=0.0
  PRIMER_WT_END_STABILITY=0.0
  PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0
  PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0
  PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0
  PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0
  PRIMER_PAIR_WT_COMPL_ANY_TH=0.0
  PRIMER_PAIR_WT_COMPL_END_TH=0.0
  PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH=0.0
  PRIMER_PAIR_WT_COMPL_ANY=0.0
  PRIMER_PAIR_WT_COMPL_END=0.0
  PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0
  PRIMER_PAIR_WT_DIFF_TM=0.0
  PRIMER_PAIR_WT_LIBRARY_MISPRIMING=0.0
  PRIMER_PAIR_WT_PR_PENALTY=1.0
  PRIMER_PAIR_WT_IO_PENALTY=0.0
  PRIMER_INTERNAL_MIN_SIZE=18
  PRIMER_INTERNAL_OPT_SIZE=20
  PRIMER_INTERNAL_MAX_SIZE=27
  PRIMER_INTERNAL_MIN_TM=57.0
  PRIMER_INTERNAL_OPT_TM=60.0
  PRIMER_INTERNAL_MAX_TM=63.0
  PRIMER_INTERNAL_MIN_GC=20.0
  PRIMER_INTERNAL_OPT_GC_PERCENT=50.0
  PRIMER_INTERNAL_MAX_GC=80.0
  PRIMER_INTERNAL_MAX_SELF_ANY_TH=47.00
  PRIMER_INTERNAL_MAX_SELF_END_TH=47.00
  PRIMER_INTERNAL_MAX_HAIRPIN_TH=47.00
  PRIMER_INTERNAL_MAX_SELF_ANY=12.00
  PRIMER_INTERNAL_MAX_SELF_END=12.00
  PRIMER_INTERNAL_MIN_QUALITY=0
  PRIMER_INTERNAL_MAX_NS_ACCEPTED=0
  PRIMER_INTERNAL_MAX_POLY_X=5
  PRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00
  PRIMER_INTERNAL_SALT_MONOVALENT=50.0
  PRIMER_INTERNAL_DNA_CONC=50.0
  PRIMER_INTERNAL_SALT_DIVALENT=1.5
  PRIMER_INTERNAL_DNTP_CONC=0.0
  PRIMER_INTERNAL_WT_SIZE_LT=1.0
  PRIMER_INTERNAL_WT_SIZE_GT=1.0
  PRIMER_INTERNAL_WT_TM_LT=1.0
  PRIMER_INTERNAL_WT_TM_GT=1.0
  PRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0
  PRIMER_INTERNAL_WT_GC_PERCENT_GT=0.0
  PRIMER_INTERNAL_WT_SELF_ANY_TH=0.0
  PRIMER_INTERNAL_WT_SELF_END_TH=0.0
  PRIMER_INTERNAL_WT_HAIRPIN_TH=0.0
  PRIMER_INTERNAL_WT_SELF_ANY=0.0
  PRIMER_INTERNAL_WT_SELF_END=0.0
  PRIMER_INTERNAL_WT_NUM_NS=0.0
  PRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0
  PRIMER_INTERNAL_WT_SEQ_QUAL=0.0
  PRIMER_INTERNAL_WT_END_QUAL=0.0
  );

  # Incorporate web params
  for my $p (@webParam){
    my($key,$value)=split(/=/,$p);
    if(defined $param{$key}){
      #logmsg "Skipping $key => $value";
    } else {
      $param{$key}=$value;
    }
  }

  # Write the parameters
  for my $key(keys(%param)){
    print $paramFh $key."=".$param{$key}."\n";
  }
  print $paramFh "=\n";
  close $paramFh;

  return $paramFile;
}

sub usage{
  "$0: run primer3 on a table of SNP positions
  Usage: $0 table.tsv in.fasta [in2.fasta...]
  --tempdir  ''   Specify if you want to make a permanent
                  intermediate directory.
  "
}
