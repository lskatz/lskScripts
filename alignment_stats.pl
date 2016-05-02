#! /usr/bin/env perl
use Bio::AlignIO;
use Data::Dumper;
use Number::Format;
use strict;
use warnings;

die "Usage: $0 alignment1.fna alignment2.clw ..." if(@ARGV<1);
exit(main(\@ARGV));

sub main{
  my($file)=@_;
  foreach my $f (@$file){
    if(!-f $f){
      print "Error: the file $f does not exist\n";
      next;
    }
    my $stats=alignmentStats($f);

    print join("\t",("File",$f))."\n";
    while(my($stat,$value)=each(%$stats)){
      print join("\t",($stat,$value));
      print "\n";
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
  my $numberFormatter=new Number::Format(-thousands_sep=>',',-decimal_point=>'.',decimal_digits=>2);

  my ($file)=@_;
  my $format=guess_format($file);
  my $return="";
  my $defaultStats={"_file"=>$file,"fileSize"=>-s $file,"no_alignments"=>0,format=>$format};
  my $stats=\%$defaultStats;
  my $in=Bio::AlignIO->new(-file=>$file,-format=>$format);
  while(my $aln=$in->next_aln){
    my $match_line=$aln->match_line;

    $$stats{no_alignments}++;
    $$stats{length}+=$aln->length;
    $$stats{is_not_flush}+=bool(!$aln->is_flush); # not flush one time will make it not flush throughout
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

  $$stats{is_flush}=!$$stats{is_not_flush};
  delete($$stats{is_not_flush});

  # format some numbers to make it easier to look at
  # TODO make an option to turn this on/off
  $$stats{fileSize}=$numberFormatter->format_bytes($$stats{fileSize},precision=>2);
  for(qw(score percentage_identity length)){
    $$stats{$_}=$numberFormatter->format_number($$stats{$_},2);
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
    #$score+=colScore(\@col,\%matrix); #LK this was producing warnings
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
  
sub bool{
  my $value=shift;
  return 1 if($value);
  return 0;
} 
