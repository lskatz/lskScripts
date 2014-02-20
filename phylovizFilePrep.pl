#!/usr/bin/env perl
# use this to prepare a SNP input file for phyloviz
# https://github.com/josephhughes/PrepPhyloViz
# author: Joseph Hughes
# edited by Lee Katz

use Bio::SeqIO;
use strict;
use warnings;
use Getopt::Long; 
use Data::Dumper;

my ($infile,$outfile,$help);
&GetOptions(
	    'in:s'      => \$infile,#input alignment in fasta format
	    'output:s'   => \$outfile,#output file for phyloviz 
      'help' => \$help,
           );
if($help || !($infile && $outfile)){
  die "Usage: $0 -i infile.fasta -o out.phyloviz\n where the first sequence is a reference sequence";
}

my $in  = Bio::SeqIO->new(-file => "$infile" , '-format' => 'fasta');
open (OUT, ">SNP$outfile")|| die "Can't open SNP$outfile\n";
open (OUT2, ">IsolateData$outfile")|| die "Can't open IsolateData$outfile\n";
my (@refbases, %bases, %mismatch, $i, %polysites);

# get bases for reference sequence
my $refseq = $in->next_seq();
@refbases=split(//,uc($refseq->seq));
$bases{$refseq->id}=join("",@refbases);

# get bases for other sequences
while ( my $seq_obj = $in->next_seq() ) {
    my $id=$seq_obj->display_id;
    $bases{$id}=uc($seq_obj->seq);
}

# find variants (SNPs)
foreach my $key(keys %bases){ 
   my @splitseq=split(//,$bases{$key});
   #print "$key $splitseq[0] $splitseq[1]\n";
   for ($i=0; $i<scalar(@splitseq); $i++){
     if ($splitseq[$i]!~/$refbases[$i]/){
       #print "mismatch at $i ref is $refbases[$i] and variant is $splitseq[$i]\n";
       $polysites{$i}++;
     }
     $mismatch{$key}{$i}=$splitseq[$i];
   }
}

my (@SNPsites, %seenSNPprofile);
my $uniqid=1;
foreach my $polysite (keys %polysites){
  #print "$polysite\n";
  push(@SNPsites,$polysite);
}
my @ascending = sort { $a <=> $b } @SNPsites;

# make unique IDs (MLST ID, or whatever) for each profile
my %uniqid;
for my $id (keys %mismatch){
  my $SNPprofile="";
  foreach my $SNP(@ascending){
    $SNPprofile.="$mismatch{$id}{$SNP}\t";
  }
  my $year="";
  $year=$1 if ($id=~m/(\d{4})/);
  $uniqid{$SNPprofile.$year}=$year."_".$uniqid++;
}

# generate output files
print OUT "ID\t";
foreach my $SNP(@ascending){
   print OUT "SNP$refbases[$SNP]".($SNP+1)."\t";
}
print OUT "\n";
print OUT2 "ID\tSeqID\n";
for my $id (keys %mismatch){
  my $SNPprofile="";
  foreach my $SNP(@ascending){
    $SNPprofile.="$mismatch{$id}{$SNP}\t";
  }
  my $year="";
  $year=$1 if ($id=~m/^(\d{4})/);
  $uniqid=$uniqid{$SNPprofile.$year} or die "Internal error: could not find unique id for $id";
  
  # link profile to genome id
  print OUT2 "$uniqid\t$id\n";
  # put the profile to the file
  if (!$seenSNPprofile{$SNPprofile}){
    print OUT "$uniqid\t$SNPprofile\n";
    $seenSNPprofile{$SNPprofile}=$uniqid;
  }
}
