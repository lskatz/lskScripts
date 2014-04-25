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
  GetOptions($settings,qw(prefix=s help acceptIds=s sitesFile=s));
  die usage() if($$settings{help});
  
  my $prefix=$$settings{prefix}||"$0.out";
  my $alnFile=shift(@ARGV) || die "ERROR: need an alignment file\n".usage();
  die "ERROR: you gave more than one alignment file: ".join(", ",@ARGV)."\n".usage() if(@ARGV>1);
  my $sitesFile=$$settings{sitesFile} || "";

  my($strainSeq)=seqToIdHash($alnFile,$settings);
  my @site=getSites($sitesFile,$strainSeq,$settings);
  #my %site; @site{@site}=@site; # generate an index of sites, shorthand
  $strainSeq=removeUninformativeSites($strainSeq,\@site,$settings);
  my($STstrains,$STseq)=sequencesToSTs($strainSeq,$settings);

  printResults($STstrains,$STseq,\@site,$prefix,$settings);

  return 0;
}

sub getSites{
  my($sitesFile,$strainSeq,$settings)=@_;
  my @site;
  my $site;
  if(-f $sitesFile){
    $site=read_file($sitesFile);
    @site=split(/\s+/,$site);
  }else{
    warn "WARNING: $sitesFile is not a positions file\n";
    my $sequence=(keys(%$strainSeq))[0];
    push(@site,$_) for(1..length($sequence));
  }
  
  return @site;
}

sub seqToIdHash{
  my($file,$settings)=@_;

  my $aln=Bio::AlignIO->new(-file=>$file)->next_aln;
  my @seq=$aln->each_seq;

  my %acceptId;
  if($$settings{acceptIds}){
    my @id=split(/\n/,`cat '$$settings{acceptIds}'`);
    $acceptId{$_}=1 for(@id); # produce warning
  }

  my %strainSeq;
  for(@seq){
    my ($genomeId,$seq)=($_->id,$_->seq);

    next if(keys(%acceptId) && !$acceptId{$genomeId});
    $strainSeq{$seq}=$genomeId;
  }

  return \%strainSeq;
}

sub sequencesToSTs{
  my($strainSeq,$settings)=@_;
  
  # Make a new ST system.
  my $STcounter=1;
  my (%STstrains,%STseq);
  while(my($genomeId,$seq)=each(%$strainSeq)){

    # assign an ST
    my $currST;
    if($currST=$STseq{$seq}){
    } else {
      $currST=$STcounter;
      $STseq{$seq}=$currST;
      $STcounter++;
    }
    push(@{ $STstrains{$currST} },$genomeId);
    $STseq{$seq}=$currST;
  }

  return (\%STstrains,\%STseq);
}

sub printResults{
  my($STstrains,$STseq,$snpSite,$prefix,$settings)=@_;

  open(ALN,">","$prefix.aln.fas") or die "Could not open aln file: $!";
  open(ALNFULL,">","$prefix.full.aln.fas") or die "Could not open full aln file: $!";
  open(ST, ">","$prefix.STs") or die "Could not open st file: $!";
  open(INFO,">","$prefix.info") or die "Could not open info file: $!";
  print INFO join("\t",qw(ST IDs))."\n";
  while(my($seq,$ST)=each(%$STseq)){
    my @genomeId=@{ $$STstrains{$ST} };
    print ALN ">$ST\n$seq\n";
    print ST join("\t",$ST,@genomeId)."\n";
    print INFO join("\t",$ST,join(",",@genomeId))."\n";

    # make a full fasta file, without STs, with duplicate sequences
    for(@genomeId){
      print ALNFULL ">$_\n$seq\n";
    }
  }
  close INFO;
  close ST;
  close ALNFULL;
  close ALN;

  open(SITES,">","$prefix.sites") or die "Could not open sites file: $!";
  for (@$snpSite){
    print SITES "$_\t" if(defined($_));
  }
  print SITES "\n";
  close SITES;

  return 1;
}

# remove uninformative sites from a STstrains hash
sub removeUninformativeSites{
  my($STstrains,$snpSites,$settings)=@_;

  # make a 2d hash of bases
  my %seq;
  my $length; # length of alignment
  while(my($seq,$ST)=each(%$STstrains)){
    my @tmp=split(//,$seq);
    $seq{$ST}=\@tmp;
    $length=@tmp;
  }

  # find sites that are nonvariant
  my %badSites;
  for(my $i=0;$i<$length;$i++){
    my $colHasVariation=0;
    my($refST,$refBasesArr)=each(%seq);
    my $refBase=$$refBasesArr[$i];
    while(my($ST,$basesArr)=each(%seq)){
      my $nt=$$basesArr[$i];
      $colHasVariation=1 if($nt ne $refBase);
      $badSites{$i}=1 if($nt=~/[nN]/);
    }
    $badSites{$i}=1 if(!$colHasVariation || $refBase=~/[nN]/);
  }

  # make a new STstrains hash, removing nonvariant sites
  my $newSeq;
  for(my $i=0;$i<$length;$i++){
    next if($badSites{$i});
    while(my($ST,$basesArr)=each(%seq)){
      $$newSeq{$ST}.=$$basesArr[$i];
    }
  }

  # remove the bad sites from the sites array
  for(my $i=0;$i<@$snpSites;$i++){
    if($badSites{$i}){
      delete($$snpSites[$i]);
    }
  }

  return $newSeq;
}

sub usage{
  local $0=fileparse($0);
  "Make an alignment suitable for PhyloViz.
  Converts an alignment into sequence types and the boiled down alignment, removing redundant ST entries.
  Usage: $0 file.aln -p prefix
    -p prefix for output files: \$p.STs.aln and \$p.STs

    optional:
    -a accept.txt a file of IDs to accept. Others will be tossed out. Default: accept all
    -s input sites file so that you can know which positions were retained
  "
}
