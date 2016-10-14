#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Bio::Perl;
use Data::Dumper;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help pos-control|positive-control first-char|first-characters=i at-least=i adapter-set=s list-adapters)) or die $!;
  die usage() if($$settings{help});

  $$settings{'first-char'}||=9999;
  $$settings{'at-least'}||=0;
  $$settings{'adapter-set'}||="skewer"; 
  $$settings{'adapter-set'}=lc($$settings{'adapter-set'});

  my(@fastq)=@ARGV;

  # echo $ADAPTER | perl -Mautodie -MBio::Perl -e '$in=Bio::SeqIO->new(-fh=>\*STDIN,-format=>"fasta"); while($seq=$in->next_seq){ my $sequence=uc($seq->seq); my %numAdapter; for my $fastq(glob("*.fastq.gz")){ open(FASTQ,"zcat $fastq | "); while(my $id=<FASTQ>){ $read=uc(<FASTQ>); $numAdapter{$sequence}++ if($read=~/$sequence/); } close FASTQ; } for(keys(%numAdapter)){print join("\t",$_,$numAdapter{$_})."\n";} last; }'
  my $adapterSeq=adapters($settings);

  for my $fastq(@fastq){
    my $numAdapter=countAdaptersInFastq($fastq,$adapterSeq);
    # Report
    while(my($adapterName,$adapterSeq)=each(%$adapterSeq)){
      my $count=$$numAdapter{$adapterName};
      next if($count < $$settings{'at-least'});
      print join("\t",$fastq,$adapterName,$adapterSeq,$count)."\n";
    }
  }

  return 0;
}

sub countAdaptersInFastq{
  my($fastq,$adapterSeq)=@_;
  my @adapterName=keys(%$adapterSeq);
  my @adapterSeq=values(%$adapterSeq);
  my $numAdapters=@adapterName;

  my %numAdapter = map{$_=>0} @adapterName;
  open(FASTQ,"zcat $fastq | ") or die "ERROR: could not read $fastq: $!";
  while(my $id=<FASTQ>){
    my $read=uc(<FASTQ>);
    for(my $i=0;$i<$numAdapters;$i++){
      if($read=~/$adapterSeq[$i]/){
        $numAdapter{$adapterName[$i]}++;
      }
    }
    # Burn two lines
    my $plus=<FASTQ>;
    my $qual=<FASTQ>;
  }
  close FASTQ;
  return \%numAdapter;
}



# Read adapters as uppercase. Return a hash
# of seqid=>sequence.
sub adapters{
  my($settings)=@_;

  my %adapterSet=(
    trimmomatic=>">NexteraPE-PrefixNX/1
AGATGTGTATAAGAGACAG
>NexteraPE-PrefixNX/2
AGATGTGTATAAGAGACAG
>NexteraPE-Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>NexteraPE-Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>NexteraPE-Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>NexteraPE-Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

>TruSeq2PEPrefixPE/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TruSeq2PEPrefixPE/2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>TruSeq2PEPCR_Primer1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TruSeq2PEPCR_Primer1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
>TruSeq2PEPCR_Primer2
CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
>TruSeq2PEPCR_Primer2_rc
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq2PEFlowCell1
TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC
>TruSeq2PEFlowCell2
TTTTTTTTTTCAAGCAGAAGACGGCATACGA

>TruSeq2_SE
AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
>TruSeq2_PE_f
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>TruSeq2_PE_r
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG

>TruSeq3PE-PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TruSeq3PE-PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

>TrueSeq3PE2-PrefixPE/1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TrueSeq3PE2-PrefixPE/2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>TrueSeq3PE2-PE1
TACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TrueSeq3PE2-PE1_rc
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>TrueSeq3PE2-PE2
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>TrueSeq3PE2-PE2_rc
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

>TruSeq3SE_IndexedAdapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>TruSeq3SE_UniversalAdapter
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
",
  # https://github.com/vsbuffalo/scythe/blob/master/illumina_adapters.fa
  skewer => ">multiplexing-forward
GATCGGAAGAGCACACGTCT
>solexa-forward
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>truseq-forward-contam
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
>truseq-reverse-contam
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
>nextera-forward-read-contam
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>nextera-reverse-read-contam
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>solexa-reverse
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
  ",
  );

  if($$settings{'list-adapters'}){
    print join("\n",keys(%adapterSet))."\n";
    die;
  }
  

  my %seq;
  my $in=Bio::SeqIO->new(-string=>$adapterSet{$$settings{'adapter-set'}},-format=>"fasta",-alphabet=>"dna");
  while(my $seq=$in->next_seq){
    $seq{$seq->id}=substr(uc($seq->seq),0,$$settings{'first-char'});
  }

  if($$settings{'pos-control'}){
    $seq{"positive-control"}="AAAA";
  }
  return \%seq;
}

sub usage{
  "$0: Detect which adapters are in your fastq file
  Usage: $0 file.fastq.gz [file2.fastq.gz...]
  --first-char   999     Takes the first X characters of the 
                         adapter sequence
  --pos-control  0       Adds in a positive control, 'AAAA'
                         to the list of primers
  --adapter-set  skewer  Which set of primers to use?
                         Options: skewer, trimmomatic

  OUTPUT
  --at-least     0       Must have at least this count to
                         be reported.
  --list-adapters        List choices of adapters sets
  "
}

