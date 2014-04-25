#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Bio::Perl;
use Bio::AlignIO;
use Getopt::Long;
use File::Basename;
use File::Slurp qw/read_file/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(prefix=s help));
  die usage() if($$settings{help});
  
  my $prefix=$$settings{prefix}||"$0.out";
  my $alnFile=shift(@ARGV) || die "ERROR: need an alignment file\n".usage();
  die "ERROR: you gave more than one alignment file: ".join(", ",@ARGV)."\n".usage() if(@ARGV>1);

  my($strainSeq)=seqToIdHash($alnFile,$settings);
  printResults($strainSeq,$prefix,$settings);

  return 0;
}

# Make a hash of sequence => [id1,id2,...]
sub seqToIdHash{
  my($alnFile,$settings)=@_;
  my $aln=Bio::AlignIO->new(-file=>$alnFile)->next_aln;
  
  my %strainSeq;
  for my $seq($aln->each_seq){
    push(@{$strainSeq{$seq->seq}},$seq->id);
  }
  return \%strainSeq;
}

sub printResults{
  my($strainSeq,$prefix,$settings)=@_;
  open(ALN,">","$prefix.aln.fas") or die "ERROR: Could not write to alignment file $prefix.aln.fas:$!";
  open(STS,">","$prefix.STs.txt") or die "ERROR: Could not write to alignment file $prefix.STs.txt:$!";
  print STS join("\t",qw(ST IDs))."\n";
  my $STcounter=0;
  while(my($sequence,$idArr)=each(%$strainSeq)){
    $STcounter++;
    print ALN ">$STcounter\n$sequence\n";
    print STS join("\t",$STcounter,join(" ",@$idArr))."\n";
  }
  close ALN; close STS;
  print "Wrote $prefix.aln.fas and $prefix.STs.txt\n";
}
sub usage{
  local $0=fileparse($0);
  "Make an alignment suitable for PhyloViz.
  Converts an alignment into sequence types and the boiled down alignment, removing redundant ST entries.
  Usage: $0 file.aln -p prefix
    -p prefix for output files: \$p.STs.aln and \$p.STs
  "
}
