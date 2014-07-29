#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Temp ('tempdir');

sub logmsg{print STDERR "@_\n";}
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help reference=s tempdir=s numcpus=i paired method=s)) or die;
  $$settings{tempdir}||=mktempdir();
  $$settings{numcpus}||=1;
  $$settings{paired}||=0;
  $$settings{reference}||="";
  $$settings{method}||="map";
  $$settings{method}=lc($$settings{method});
  
  die usage() if($$settings{help});

  my ($fastq1,$fastq2)=writeFastq($settings);

  my ($reads1,$reads2)=({},{}); # hashes to convert to fastq, to write to stdout
  if($$settings{method} eq 'map'){
    ($reads1,$reads2)=sortByMapping($$settings{reference},$fastq1,$fastq2,$settings);
  } elsif($$settings{method} eq 'binbygc' || $$settings{method} eq 'binbyalpha'){
    ($reads1,$reads2)=sortByBinning($fastq1,$fastq2,$settings);
  } else{
    die "Internal error: could not figure out what kind of sorting you wanted (user specified '$$settings{method}')";
  }


  return 0;
}

sub sortByMapping{
  my($reference,$fastq1,$fastq2,$settings)=@_;
  if(!$reference){
    ...;
    logmsg "A reference was not specified. Assembling the reads.";
    $reference=assembleReads([$fastq1,$fastq2],$settings);
  } 
  die "ERROR could not find reference genome at $reference\n".usage() if(!-f $reference);

  my $bam=mapReads($reference,$fastq1,$fastq2,$settings);
  my $reads1=readFastq($fastq1,$settings);
  my $reads2=readFastq($fastq2,$settings);
  my $numReads=printReadsFromBam($bam,$reads1,$reads2,$settings);
  logmsg "$numReads reads or pairs printed to stdout";
}

sub writeFastq{
  my($settings)=@_;

  my $is_paired=$$settings{paired};

  # write the fastq
  my $newfastq1="$$settings{tempdir}/input1.fastq";
  my $newfastq2="$$settings{tempdir}/input2.fastq";
  logmsg "Reading STDIN to create a fastq file";
  open(FASTQ1,">",$newfastq1) or die "ERROR: could not open temporary fastq file at $newfastq1: $!";
  open(FASTQ2,">",$newfastq2) or die "ERROR: could not open temporary fastq file at $newfastq2: $!" if($is_paired);
  system("echo '' > $newfastq2") if(!$is_paired); die if $?;
  my $i=0;
  while(<>){
    my $mod=$i++ % 8;
    if($mod<4 || !$is_paired){
      print FASTQ1 $_;
    } else {
      print FASTQ2 $_;
    }
  }
  close FASTQ1;
  close FASTQ2 if($is_paired);
  return ($newfastq1,$newfastq2) if wantarray;
  return $newfastq1;
}

sub assembleReads{
  my($fastq,$settings)=@_;
  my $assembly="$$settings{tempdir}/assembly.fasta";
  system("run_assembly_illumina.pl --fast '$fastq' -o $assembly 1>&2");
  die if $?;
  return $assembly;
}

sub mapReads{
  my($reference,$fastq1,$fastq2,$settings)=@_;
  my $sam="$$settings{tempdir}/smalt.sam";
  my $bam="$$settings{tempdir}/smalt.bam";
  my $sortedPrefix="$$settings{tempdir}/smalt.sorted";
  my $sorted="$sortedPrefix.bam";

  system("smalt index '$reference' '$reference'"); die if $?;
  my $command="smalt map -f samsoft -n $$settings{numcpus} '$reference' '$fastq1' ";
     $command.="'$fastq2' " if($$settings{paired});
     $command.="> $sam";
  system($command); die if $?;
  logmsg "Done mapping. Transforming with Samtools";
  system("samtools view -bS '$sam' -T '$reference' > '$bam'"); die if $?;
  system("samtools sort '$bam' '$sortedPrefix'"); die if $?;
  system("samtools index '$sorted'"); die if $?;
  return $sorted;
}

sub mktempdir(;$) {
  my ($settings) = @_;
  my $tempdir_path = File::Spec->join(File::Spec->tmpdir(), (split("::",(caller(1))[3]))[1].".$$.XXXXX");
  my $tempdir = tempdir($tempdir_path, CLEANUP => !($$settings{keep}));
  logmsg "Created tempdir $tempdir";
  return $tempdir;
}

sub printReadsFromBam{
  my($bam,$reads1,$reads2,$settings)=@_;
  my $i=0;
  if($$settings{paired}){
    open(BAM,"samtools view -f 0x40 '$bam' | ") or die "ERROR: Could not open $bam: $!";
  } else {
    open(BAM,"samtools view '$bam' | ") or die "ERROR: Could not open $bam: $!";
  }

  my %seen;
  while(<BAM>){
    chomp;
    my @F=split /\t/;
    if(!$$reads1{$F[0]}){
      print Dumper [$reads1,\@F];die;
    }
    die "ERROR: I have already seen ID $F[0] and so this might be a shuffled paired end file" if($seen{$F[0]}++);
    print $$reads1{$F[0]};
    print $$reads2{$F[0]} if($$settings{paired});
    $i++;
  }
  return $i;
}

sub sortByBinning{
  my($fastq1,$fastq2,$settings)=@_;
  my $unsorted1=readFastqIntoHash($fastq1,$settings);
  my $unsorted2=readFastqIntoHash($fastq2,$settings);

  my @sorted1;
  if($$settings{method} eq 'binbyalpha'){
    @sorted1=sort{$$a{seq} cmp $$b{seq} } values(%$unsorted1);
  } elsif ($$settings{method} eq 'binbygc'){
    @sorted1=sort{$$a{gc} <=> $$b{gc}} values(%$unsorted1);
  } else {
    die "ERROR: I could not understand method ($$settings{method})";
  }

  # print out the results
  my $numReads=@sorted1;
  for(my $i=0;$i<$numReads;$i++){
    print "$sorted1[$i]{defline}\n$sorted1[$i]{seq}\n+\n$sorted1[$i]{qual}\n";
    if($$settings{paired}){
      my $key=$sorted1[$i]{key};
      print "$$unsorted2{$key}{defline}\n$$unsorted2{$key}{seq}\n+\n$$unsorted2{$key}{qual}\n" if($$settings{paired});
    }
  }
}

sub readFastq{
  my($file,$settings)=@_;
  logmsg "Reading $file into memory";
  my $seqs;
  open(FASTQ,$file) or die "ERROR: cannot read $file: $!";
  while(my $id=<FASTQ>){
    chomp $id;
    my ($key,$desc)=split(/\s+/,$id);
    next if(!$key);
    $key=~s/^\@//;
    $$seqs{$key}.="$id\n".<FASTQ>; # sequence and identifier
    <FASTQ>; # burn the + line
    $$seqs{$key}.="+\n"; # but add in my own
    $$seqs{$key}.=<FASTQ>; # qual
  }
  close FASTQ;
  return $seqs;
}
sub readFastqIntoHash{
  my($file,$settings)=@_;
  logmsg "Reading $file into memory";
  my $seqs;
  open(FASTQ,$file) or die "ERROR: cannot read $file: $!";
  while(my $defline=<FASTQ>){
    my ($key,@desc)=split(/\s+/,$defline);
    die "ERROR: no key for $defline" if(!$key);
    $key=~s/^\@//;
    $$seqs{$key}{defline}=$defline; # whole defline
    $$seqs{$key}{key}=$key;         # text before whitespace
    $$seqs{$key}{desc}=join(" ",@desc);
    $$seqs{$key}{seq}=<FASTQ>;
    <FASTQ>; # burn the + line
    $$seqs{$key}{qual}=<FASTQ>;
    chomp($$seqs{$key}{qual},$$seqs{$key}{defline},$$seqs{$key}{seq},$$seqs{$key}{desc});

    # some derivative values that will need to be calculated
    $$seqs{$key}{length}=length($$seqs{$key}{seq});
    $$seqs{$key}{gc}=($$seqs{$key}{seq}=~s/([GC])/$1/gi)/$$seqs{$key}{length};
  }
  close FASTQ;
  return $seqs;
}

sub readFastqDefline{
  my($file,$settings)=@_;
  my $defline;
  open(FASTQ,$file) or die "ERROR: cannot read $file: $!";
  while(my $id=<FASTQ>){
    chomp $id;
    $id=~s/^\@//;
    my($identifer,@desc)=split(/\s+/,$id);
    my $desc=join(" ",@desc);
    <FASTQ> for(1..3); # burn these lines
    $$defline{$identifer}=$desc;
  }
  close FASTQ;
  return $defline;
}

sub usage{
  "Allows for better compression of reads by sorting them. Uses smalt under the hood. CGP fast assembly for an assembler.
  Usage: $0 [-r reference.fasta] < reads.fastq | gzip -c > sorted.fastq.gz
  -r reference.fasta If not supplied, then the reads will be assembled to create one.
  -m methodology  Values: map, binByAlpha, binByGC  Default: map
    map: map reads onto a reference genome given by -r
    binByAlpha: alphabetize the reads
    binByGC: bin the reads according to their nucleotide composition
  -t tempdir/ Default: one will be created for you
  --numcpus How many cpus to use when mapping. Default: 1
  -p to describe paired-end
  "
}
