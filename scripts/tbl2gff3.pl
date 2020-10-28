#!/usr/bin/env perl

use Getopt::Long qw/GetOptions/;
use Data::Dumper qw/Dumper/;
use strict;
use warnings;

#####################################
#	Program: Tbl2KnownGene
#	Authors: Yongsheng Bai
#	Modified by: Lee Katz
#   Copyright 2014
#	Description: Tbl2KnownGene is a command-line program that parses the contents of a NCBI .tbl annotation table and produces a UCSC Known Genes annotation table. Tbl2KnownGene is written in Perl.
#	Command: perl Tbl2KnownGene.pl NCBI_Chr1.tbl NCBI_Chr2.tbl NCBI_Chr3.tbl NCBI_Chr4.tbl NCBI_Chr5.tbl
#####################################

die "TODO - I have only just begun working on CDS";

sub logmsg{print STDERR "$0: @_\n";}
exit main();

sub main{

  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  usage() if($$settings{help});

  print "##gff-version 3\n";
  
  # http://uswest.ensembl.org/info/website/upload/gff3.html
  print "#".join("\t", qw(seqid source type start end score strand phase attributes))."\n";
  print "\n";

  my %locusTag=();
  my $prot_counter = 0;
  my $protein_id;
  my $transcript_id;
  my $flag;
  my $TRUE = 1;
  my $FALSE = -1;

  # LK adding global variables after adding 'strict'
  my $zone = "";
  my $genestart = "<1";
  my $geneend   = "inf";
  my $exon_start = "<1";
  my $exon_stop  = "inf";
  my $exon_stop_reverse  = "inf";
  my $exon_start_reverse = "inf";
  my $mrna_start = "<1";
  my $mrna_end   = "inf";
  my $cds_start = "<1";
  my $cds_end   = "inf";
  my $strand    = "";
  my $locusTag  = "UNKNOWNTAG";
  my $chr       = "CHROMOSOME";

  my $exon_count=0;
 
  while (my $line = <>) {
    chomp($line);
    my @linearray=split(/\t/,$line);
    push(@linearray, ("") x 20); # ensure this is not undefined

    if(substr($line,0,1) eq '>'){
      if($line =~ /feature/i){
        @linearray = split(/\s+/, $line);
        $chr = $linearray[1];
      }
    }

    ############ gene set
    if ($linearray[2] eq 'gene') {
      my $zone='GENE';
      $genestart=$linearray[0];
      $geneend=$linearray[1];
      if ($genestart > $geneend) {
        $strand='-';
      }
      else {
        $strand='+';
      }
    }
    if ( ($linearray[2] eq '') && ($zone eq 'GENE') ) {
      if ($linearray[3] eq 'locus_tag') {
        $locusTag = $linearray[4];
      }
      if ($linearray[3] eq 'gene_syn') {
        $locusTag = $linearray[4];
      }
    }

    #########   CDS Set
    if ($linearray[2] eq 'CDS') {
      $prot_counter=0;
      $cds_start=0;
      $cds_end=0;
      $zone='CDS';
      if($strand eq '+')
      {
        $cds_start = $linearray[0];
        $cds_end = $linearray[1];
        print  join("\t", $transcript_id, lc($chr),
                          $strand, $mrna_start, $mrna_end,
                          $cds_start, $cds_end,
                          $exon_count, $exon_start_reverse,
                          $exon_stop_reverse, $protein_id, $locusTag
                   );
        print "\n";
      }
      else#negative strand
      {
        $cds_start = $linearray[1];
        $cds_end = $linearray[0];
      }
    }
    
    #####  for rest of CDS lines
    if (($linearray[2] eq '') && ($zone eq 'CDS') && ($linearray[3] eq '') ) {
      if($strand eq '+')
      {
        $cds_end = $linearray[1];
      }
      else#negative strand
      {
        $cds_start = $linearray[1];
      }
    }
    if ($linearray[2] eq 'mRNA') {
      $exon_count=0;
      $mrna_start=0;
      $mrna_end=0;
      $exon_start="";
      $exon_stop="";
      $exon_start_reverse="";
      $exon_stop_reverse="";
      $zone='mRNA';
      $flag = $TRUE;
      if($strand eq '+')
      {
        $mrna_start = $linearray[0];
        $mrna_end = $linearray[1];
        $exon_start .= "$linearray[0]".",";
        $exon_stop .= "$linearray[1]".",";
      }
      else#negative strand
      {
        $mrna_start = $linearray[1];
        $mrna_end = $linearray[0];
        $exon_start .= "$linearray[1]".",";
        $exon_stop .= "$linearray[0]".",";
      }
      $exon_count++;
    }
    
    ##### for rest of mRNA lines
    if (($linearray[2] eq '') && ($zone eq 'mRNA') && ($linearray[3] eq '') ) {
      if($strand eq '+')
      {
        $mrna_end = $linearray[1];
        $exon_start .= "$linearray[0]".",";
        $exon_stop .= "$linearray[1]".",";
      }
      else#negative strand
      {
        $mrna_start = $linearray[1];
        $exon_start .= "$linearray[1]".",";
        $exon_stop .= "$linearray[0]".",";
      }
      $exon_count++;
    }

    if ( ($linearray[2] eq '') && ($linearray[3] ne '') ){
      if ($linearray[3] eq 'protein_id') {
        $protein_id=$linearray[4];
      }
      if ($linearray[3] eq 'transcript_id') {
        $transcript_id=$linearray[4];
        $prot_counter++;
      }
    }
    if(($prot_counter == 2)&&($flag == $TRUE))
    {
      my @exonStartArray = split(",",$exon_start);
      my @exonStopArray = split(/\,/,$exon_stop);
      if($strand eq '-')
      {
        for(my $i=($#exonStartArray); $i >=0; $i--)
        {
          $exon_start_reverse .= "$exonStartArray[$i]".",";
          $exon_stop_reverse .= "$exonStopArray[$i]".",";
        }
        print  "$transcript_id","\t",lc($chr),"\t",$strand,"\t",$mrna_start,"\t",$mrna_end,"\t",$cds_start,"\t",$cds_end,"\t",$exon_count,"\t",$exon_start_reverse,"\t",$exon_stop_reverse,"\t", $protein_id, "\t", $locusTag,"\n";
      }
      else# + strand
      {
        print  "$transcript_id","\t",lc($chr),"\t",$strand,"\t",$mrna_start,"\t",$mrna_end,"\t",$cds_start,"\t",$cds_end,"\t",$exon_count,"\t",$exon_start,"\t",$exon_stop,"\t", $protein_id, "\t", $locusTag,"\n";
      }
      $flag = $FALSE;
    }
  }

  return 0;
}

sub usage{
  print "$0 < in.tbl > out.gff\n";
  exit 0;
}

