#!/usr/bin/env perl
# Compares gene predictions - presence/absence and also start position congruence
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::Perl;
use Bio::Tools::GuessSeqFormat;
use File::Basename qw/fileparse basename dirname/;

sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help includeRef));

  my($ref,@query)=@ARGV;
  die usage() if(!$ref || !$query[0]);
  for($ref,@query){
    die "ERROR: I can't find $_\n".usage() if(!-e $_);
  }

  # Print the header and then get onto the comparisons
  print join("\t",qw(query genesInCommon commonStartSites exactSame))."\n";
  unshift(@query,$ref) if($$settings{includeRef});
  for my $query(@query){
    compareGenePredictions($ref,$query,$settings);
  }
}

sub compareGenePredictions{
  my($ref,$query,$settings)=@_;
  
  my $refPredictions=genePredictions($ref,$settings);
  my $queryPredictions=genePredictions($query,$settings);

  # Small sanity check on whether there are matching contigs
  my %refSeqnames;
  @refSeqnames{(keys(%$refPredictions))}=(1) x scalar((keys(%$refPredictions)));
  for(keys(%$queryPredictions)){
    if(!$refSeqnames{$_}){
      logmsg "SKIPPING: could not find seqname $_ in $ref, but it exists in $query!";
      return 0;
    }
  }

  # Compare presence/absence and starts
  my $presenceCount=0;
  my $startConcordance=0;
  my $totalReferenceGenes=0;
  while(my($seqname,$strands)=each(%$refPredictions)){
    while(my($strand,$CDS)=each(%$strands)){
      while(my($end,$refStart)=each(%$CDS)){
        $totalReferenceGenes++;

        # Is the gene present?
        my $queryStartSite=$$queryPredictions{$seqname}{$strand}{$end};
        next if(!defined($queryStartSite));
        $presenceCount++;

        # Does the gene have the same start site?
        next if($queryStartSite != $refStart);
        $startConcordance++;
      }
    }
  }

  # Turn these counts into percentages
  my $percentPresent=sprintf("%0.2f",($presenceCount/$totalReferenceGenes*100));
  my $percentStartConcordancy=sprintf("%0.2f",($startConcordance/$presenceCount*100));
  my $percentStartOverall=sprintf("%0.2f",($startConcordance/$totalReferenceGenes*100));

  print join("\t",$query,"$percentPresent%","$percentStartConcordancy%","$percentStartOverall%")."\n";
}

# Return a complex hash of gene predictions. Each prediction is based on the contig and stop of a gene.
#    $hash->{contig}->{strand}->{stop} = start
sub genePredictions{
  my($file,$settings)=@_;

  my $format=Bio::Tools::GuessSeqFormat->new(-file=>$file)->guess;
  if(!defined $format){
    die "ERROR: I could not guess the format of $file";
  }
  
  my $predictions={};
  if($format=~/genbank|embl/i){
    $predictions=genePredictionsFromRichseq($file,$format,$settings);
  } elsif($format=~/gff/i){
    $predictions=genePredictionsFromGff($file,$format,$settings);
  } else {
    die "I haven't yet figured out a way to get predictions out of format $format from file $file";
  }

  return $predictions;
}

sub genePredictionsFromRichseq{
  my($file,$format,$settings)=@_;
  
  my $predictions={};
  my $in=Bio::SeqIO->new(-file=>$file,-format=>$format);
  while(my $seq=$in->next_seq){
    for my $feat($seq->get_SeqFeatures){
      next if($feat->primary_tag ne "CDS");
      $$predictions{$seq->id}{$feat->strand}{$feat->end}=$feat->start;
    }
  }
  return $predictions;
}

sub genePredictionsFromGff{
  my($file,$format,$settings)=@_;

  my $predictions={};
  my $in;

  open(GFF,$file) or die "ERROR: could not open GFF $file: $!";
  while(<GFF>){
    s/^\s+|\s+$//g; # trim
    next if(/^$/);  # empty lines
    next if(/^#/);  # comments
    my($seqname, $source, $type, $start, $end, $score, $strand, $phase,$attributes)=split(/\t/);
    next if($type ne "CDS");
    $seqname=~s/\s+.*$//; # seqname should just be the first "word"
    $strand=1 if($strand=~/\+/);
    $strand=-1 if($strand=~/\-/);
    $$predictions{$seqname}{$strand}{$end}=$start;
  }
  return $predictions;
}

sub usage{
  local $0=basename $0;
  "$0: compares the presence/absence of genes across gene predictions and their start sites.
  Usage: $0 ref.gbk query.gbk [query2.gbk ...]
  --includeRef  Include the reference genome against itself to have a sanity '100%' line

  NOTES:
  Only CDS features are compared, and the accepted formats are GenBank, EMBL, and GFF.
  Contig names must be the same across comparisons.
  "
}
