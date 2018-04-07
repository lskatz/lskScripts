#! /usr/bin/env perl
use Bio::AlignIO;
use Data::Dumper;
use strict;
use warnings;
use Getopt::Long;

# globals
my @header=qw(File no_sequences no_alignments length percentage_identity format snps conserved_sites score);

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help));

  my @file=@ARGV;
  die usage() if(!@file || $$settings{help});

  print join("\t",@header)."\n";
  foreach my $f (@file){
    if(!-f $f){
      print "Error: the file $f does not exist\n";
      next;
    }
    my $stats=alignmentStats($f);
    for my $h(@header){
      die "Internal error: $h was not calculated" if(!defined($$stats{$h}));
      print $$stats{$h}."\t";
    }
    print "\n";
  }
  return 0;
}

sub guess_format{
  my $filename=shift;
  for my $format (qw(xmfa fasta pfam selex stockholm prodom clustalw msf mase bl2seq nexus phylip)){
    eval{
      my $in=Bio::AlignIO->new(-file=>$filename,-format=>$format);
      my $aln=$in->next_aln;
      $aln->length;
    };
    if($@){
      my $error=$@;
      next;
    }
    return $format;
  }
  die "Could not guess the format of $filename\n";
}

sub alignmentStats{
  my ($file)=@_;
  my $format=guess_format($file);
  my $return="";
  my $defaultStats={File=>$file, "_file"=>$file,"fileSize"=>-s $file,"no_alignments"=>0,format=>$format};
  my $stats=\%$defaultStats;
  my $in=Bio::AlignIO->new(-file=>$file,-format=>$format);
  while(my $aln=$in->next_aln){
    my $match_line=$aln->match_line;

    $$stats{no_alignments}++;
    $$stats{length}+=$aln->length;

    $$stats{no_sequences}+=$aln->num_sequences;
    $$stats{percentage_identity}+=$aln->percentage_identity;
    $$stats{score}+=score($aln);

    my $conserved_sites=($match_line=~s/\*//g);
    $$stats{conserved_sites}+=$conserved_sites;
  }
  $$stats{snps}=$$stats{length}-$$stats{conserved_sites};
  return $defaultStats if($$stats{no_alignments} == 0);

  # average out some stats.  
  # TODO This should be normalized better by weighting bigger alignments.
  foreach (qw(percentage_identity)){
    $$stats{$_}/=$$stats{no_alignments};
  }

  # round some numbers
  for (qw(percentage_identity)){
    $$stats{$_}=sprintf("%0.2f",$$stats{$_});
  }

  return $stats;
}

# assume DNA for right now
# TODO include a matrix for proteins
sub score{
  my($aln)=@_;
  my(
    $score,      # MSA score
    @sequence,   # the sequence of each Seq
  );
  $score=0;
  my %matrix=dnaMatrix();
  my $alnLength=$aln->length;
  my @seq=$aln->each_seq;
  push(@sequence,$_->seq) for(@seq);
  for(my $i=0;$i<$alnLength;$i++){
    my @col=();
    for (@sequence){
      my $nt=substr($_,$i,1);
      push(@col,$nt);
    }
    $score+=colScore(\@col,\%matrix);
  }

  return $score;
}

sub colScore{
  my($col,$matrix)=@_;
  my $s=0;
  my $numLetters=@$col;
  for(my $i=0;$i<$numLetters;$i++){
    my $letter1=$$col[$i];
    for(my $j=$i+1;$j<$numLetters;$j++){
      my $letter2=$$col[$j];
      $s+=$$matrix{$letter1}{$letter2};
    }
  }
  return $s;
}

sub dnaMatrix{
  my $m=2;    # match
  my $mm=-1;  # mismatch
  my $g=-1;   # gap
  my %matrix=(
    A=>{A=>$m,T=>$mm,C=>$mm,G=>$mm,N=>$mm,'-'=>$g}, 
    T=>{A=>$mm,T=>$m,C=>$mm,G=>$mm,N=>$mm,'-'=>$g},
    C=>{A=>$mm,T=>$mm,C=>$m,G=>$mm,N=>$mm,'-'=>$g},
    G=>{A=>$mm,T=>$mm,C=>$mm,G=>$m,N=>$mm,'-'=>$g},
    N=>{A=>$mm,T=>$mm,C=>$mm,G=>$mm,N=>$mm,'-'=>$g}, # all mm for N
    '-'=>{A=>$g,T=>$g,C=>$g,G=>$g,N=>$g,'-'=>$g},
  );
  return %matrix;
} 

sub usage{
  "Usage: $0 alignment1.fna alignment2.clw ...
  "
}
