#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;

use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Align::DNAStatistics;
use Bio::Tree::DistanceFactory;
use Bio::TreeIO;

local $0=basename $0;
sub logmsg{print STDERR "@_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings, qw(help scheme=s)) or die $!;
  die usage() if(!$$settings{scheme} || !@ARGV);

  my($mlstProfile)=@ARGV;

  logmsg "Generating the fasta alignment of MLST alleles";
  my $fasta=makeFasta($mlstProfile, $$settings{scheme}, $settings);
  logmsg "Generating the tree";
  my $tree=makeTree($fasta, $settings); # TODO make a tree with bioperl
  
  # Print to stdout
  my $treeout = Bio::TreeIO->new(-format=>"newick");
  $treeout->write_tree($tree);
  print "\n";

  return 0;
}

sub makeFasta{
  my($mlstProfile, $scheme, $settings)=@_;

  # Get the allele database
  my $allele=alleleDatabase($scheme,$settings);
  #die Dumper $allele;
  # Locus names are those that don't start with a capital letter
  my @locus = grep {!/^[A-Z]/} keys(%$allele);

  # Read the profiles from mlst
  open(my $fh, $mlstProfile) or die "ERROR: could not open $mlstProfile: $!";
  my $header=<$fh>; chomp $header;
  my @header=split(/\t/,$header);
  my $fasta="";
  while(<$fh>){
    chomp;
    my %F;
    @F{@header} = split /\t/;
    
    # go from low-confidence to high confidence
    for(@locus){
      $F{$_}=~s/^~//;  # remove tildes
      $F{$_}=~s/\?$//; # remove question mark
    }
    
    my $fastaEntry=">".basename($F{FILE})."\n";
    for my $locus(@locus){
      $fastaEntry.=$$allele{$locus}{$F{$locus}}
    }
    #$fastaEntry=~s/(.{60})/$1\n/g; # put in newlines every 60 char

    $fasta.=$fastaEntry."\n";
  }
  close $fh;

  return $fasta;
}

sub alleleDatabase{
  my($scheme,$settings)=@_;
  # allele database reading
  my %allele;
  my @profileFiles=glob("$scheme/*");
  if(!@profileFiles){
    die "ERROR: could not find any files under $scheme";
  }

  for my $fna(@profileFiles){
    eval{
      my $in=Bio::SeqIO->new(-format=>"fasta", -file=>$fna);
      $in->next_seq;
    };
    if($@){
      logmsg "Skipping $fna";
      next;
    }
    
    my $in=Bio::SeqIO->new(-format=>"fasta", -file=>$fna);
    while(my $seq=$in->next_seq){
      my($locus, $allele);
      if($seq->id =~ /(.+)[_-](\d+)/){
        $locus=$1;
        $allele=$2;
      } else { 
        next;
      }
      $allele{$locus}{$allele}=$seq->seq;
    }

    # Make the blank allele
  }

  my @locus = grep {!/^[A-Z]/} keys(%allele);
  for(@locus){
    # the default allele will be the first or second allele, for laziness sake
    my $defaultAllele=$allele{$_}{1} || $allele{$_}{2};
    $allele{$_}{'-'} = '-' x length($defaultAllele);
  }

  return \%allele;
}

sub makeTree{
  my($fasta,$settings)=@_;

  #my $fastaFh = IO::String->new($fasta);
  #open($fastaFh, "<", $fasta) or die "ERROR: reading fasta string: $!";

  my $aln = Bio::AlignIO->new(-string=>$fasta, -format=>"fasta")->next_aln;
  my $dfactory = Bio::Tree::DistanceFactory->new(-method=>"NJ");
  my $stats = Bio::Align::DNAStatistics->new;
  
  my $mat = $stats->distance(-method=>"Kimura", -align=>$aln);
  my $tree = $dfactory->make_tree($mat);

  return $tree;
}

sub usage{
  "$0: creates a tree using mlst profiles directly from the tseemann mlst software.
  Usage: $0 --scheme scheme mlst.tsv > tree.dnd
    Where scheme is a folder with the MLST scheme.
    Run `mlst --help` for more information about the datadir parameter
  "
}
