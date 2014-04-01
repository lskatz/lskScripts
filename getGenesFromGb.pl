#!/usr/bin/env perl
# get the genes out of a genbank file
# spit out a fasta

use Bio::Perl;
use Data::Dumper;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $settings={};
die usage() if(@ARGV<1);
GetOptions($settings,qw(protein infile=s locusTags geneTags trust product)) or die $!;
my $translate=$$settings{protein} || 0;
my $useLocusTags=$$settings{locusTags} || 0;
my $useGeneTags=$$settings{geneTags} || 0;
my $inFilename=$$settings{infile} or die usage();

my $seqin=Bio::SeqIO->new(-file=>$inFilename,-format=>"genbank");
my $seqout=Bio::SeqIO->new(-format=>"fasta");
my ($strain,$dir,$ext)=fileparse($inFilename,(qw(.gb .gbk)));
$strain=~s/_/-/g; # if I am gluing together by _, then I can't use it in strain names
my %seen; # keep track of genes that have been seen
# for each contig
my $numNonsenseMutants=0;
while(my $contig=$seqin->next_seq){
  my $contigName=$contig->id;
  #$contigName=~s/_.+$//;
  $contigName=~s/_/-/g; # if I am gluing together by _, then I can't use it in contig IDs
  my @feat= grep {$_->primary_tag=~/cds|gene/i} $contig->get_SeqFeatures;

  my @sortedfeat=sort{
    my $has_translation=$a->has_tag('translation');
    my $translation=($has_translation)?($a->get_tag_values('translation'))[0] : "";
    return -1 if($a->primary_tag =~/cds/i);
    return -1 if($has_translation);
    return 0;
  } @feat;
  
  foreach my $feat ( @sortedfeat ){
    my $start=$feat->start;
    my $stop=$feat->end;
    my $seq=$feat->seq;
    my $strand=$feat->location->strand;
    my (@newId);
    my $locus="";
    my $gene="";
    my $product="";
    if($feat->has_tag('locus_tag')){
      $locus=join("_",$feat->get_tag_values('locus_tag'));
    }
    if($feat->has_tag('gene')){
      $gene=join("_",$feat->get_tag_values('gene'));
    }
    if($feat->has_tag('product')){
      $product=join("_",$feat->get_tag_values('product'));
    }

    @newId=(join("_",($strain,$contigName,$start,$stop)));
    
    my $newId=join("_",@newId);
    # remove any empties in the @newId
    $newId=~s/^_+|_+$//g;
    $newId=$locus if($locus && $useLocusTags); # keep the locus name if it exists
    $newId=$gene if($gene && $useGeneTags);
    $newId=$product if($product && $$settings{product});

    $newId=~s/\s+/_/g;
    $seq->id("lcl|$newId");

    next if($seen{$seq->id}++); # don't write this seq if it's been seen.  Increment in the same step

    if($$settings{trust} && $translate && !$feat->has_tag("translation")){
      warn "You want to trust this sequence's given translation, but it doesn't exist. Skipping.\n";
      warn "  sequence ".$seq->id."\n";
      warn "  But maybe it shows up again in the feature table, so it might still be added\n";
      $seen{$seq->id}--;
      next;
    }

    if(!$$settings{trust} && $seq->translate->seq=~/\*\w/){
      warn "There is an internal stop codon for ".$seq->id.". Attempting revcom.\n";
      $seq=$seq->revcom;
      if($seq->translate->seq=~/\*\w/){
        warn "  ". $seq->id." has an internal stop codon still. Reverting.\n";
        $seq=$seq->revcom;
        $numNonsenseMutants++;
      } else { warn "  Success!\n"; }
    }
    $seq=$seq->translate if($translate);
    # set the sequence to the feature tag if it exists and if "trust" is set
    if($$settings{trust} && $translate && $feat->has_tag('translation')){
      warn "Trusting and setting the seq to the given translation\n";
      my $translation=($feat->get_tag_values('translation'))[0];
      $seq->seq($translation);
    }

    $seqout->write_seq($seq);
  }
}

warn "There were $numNonsenseMutants genes with internal stop codons.\n" if($numNonsenseMutants);

sub usage{
  my $thisFile=fileparse($0);
  "Extracts genes from a rich seq file and prints the genes to stdout
  Usage: $thisFile -i infile
  -i infile
    genbank or embl
  --protein
    if you would like to extract protein
  -l
    to use locus tags instead of genome_contig_start_stop format
  -g
    to use gene tags instead of anything else (overrides -l option)
  -t
    to trust a given translation, if in the features (if using -p) 
  --product
    to trust the product name
  ";
}
