#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Temp qw/tempdir/;
use File::Basename qw/basename/;
use threads;
use Thread::Queue;

local $0=basename $0;
my @fastqExt=qw(.fastq.gz .fastq .fq.gz .fq);

sub logmsg{ print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help numcpus=i genomesize=i mindepth=i)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{genomesize}||=5000000;
  $$settings{mindepth}||=2;
  $$settings{tempdir}||=tempdir("MASHTREE.XXXXXX",CLEANUP=>1,TMPDIR=>1);
  logmsg "Temporary directory will be $$settings{tempdir}";

  die usage() if($$settings{help});

  my @reads=@ARGV;
  die usage() if(@reads < 2);

  # Check for prereq executables.
  for my $exe(qw(mash fneighbor)){
    system("$exe -h > /dev/null 2>&1");
    die "ERROR: could not find $exe in your PATH" if $?;
  }

  my $mshDir=sketchAll(\@reads,$settings);

  my $distances=mashDistance($mshDir,$settings);

  my $treeContent = distancesToTree($distances,$settings);

  print $treeContent;

  return 0;
}

# Run mash sketch on everything, multithreaded.
sub sketchAll{
  my($reads,$settings)=@_;

  my $sketchDir="$$settings{tempdir}/msh";
  mkdir $sketchDir;

  my $readsQ=Thread::Queue->new(@$reads);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&mashSketch,$sketchDir,$readsQ,$settings);
  }
  
  $readsQ->enqueue(undef) for(@thr);
  $_->join for(@thr);

  return $sketchDir;
}

# Individual mash sketch
sub mashSketch{
  my($sketchDir,$Q,$settings)=@_;
  while(defined(my $fastq=$Q->dequeue)){
    logmsg "Sketching $fastq";
    my $outPrefix="$sketchDir/".basename($fastq);
    if(-e "$outPrefix.msh"){
      logmsg "WARNING: ".basename($fastq)." was already mashed. You need unique filenames for this script. This file will be skipped: $fastq";
      next;
    }
    system("mash sketch -k 21 -s 10000 -m $$settings{mindepth} -c 10 -g $$settings{genomesize} -o $outPrefix $fastq > /dev/null 2>&1");
    die if $?;
  }
}

# Parallelized mash distance
sub mashDistance{
  my($indir,$settings)=@_;

  my @msh=glob("$indir/*.msh");

  my $outdir="$$settings{tempdir}/dist";
  mkdir($outdir);

  my $mshQueue=Thread::Queue->new(@msh);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&mashDist,$outdir,$mshQueue,\@msh,$settings);
  }

  $mshQueue->enqueue(undef) for(@thr);
  $_->join for(@thr);

  system("cat $outdir/*.tsv > $$settings{tempdir}/distances.tsv");
  
  return "$$settings{tempdir}/distances.tsv";
}

# Individual mash distance
sub mashDist{
  my($outdir,$mshQueue,$allMsh,$settings)=@_;
  my $allMshOpt=join(" ",@$allMsh);
  while(defined(my $msh=$mshQueue->dequeue)){
    my $outfile="$outdir/".basename($msh).".tsv";
    logmsg "Distances for $msh";
    system("mash dist -t $msh $allMshOpt > $outfile");
    die if $?;
  }
}

# 1. Read the mash distances
# 2. Create a phylip file
# 3. Run emboss's fneighbor
sub distancesToTree{
  my($distances,$settings)=@_;

  open(MASHDIST,"<",$distances) or die "ERROR: could not open $distances for reading: $!";

  my $id="UNKNOWN"; # Default ID in case anything goes wrong
  my %m; #matrix for distances
  while(<MASHDIST>){
    chomp;
    if(/^#query\s+(.+)/){
      $id=basename($1,@fastqExt);
    } else {
      my @F=split(/\t/,$_);
      $F[0]=basename($F[0],@fastqExt);
      $m{$id}{$F[0]}=sprintf("%0.6f",$F[1]);
    }
  }
  close MASHDIST;

  # Create the phylip file.
  # Make the text first so that we can edit it a bit.
  my %seenTruncName;
  my $phylipText="";
  my @genome=sort{$a cmp $b} keys(%m);
  for(my $i=0;$i<@genome;$i++){ 
    my $name=substr($genome[$i],0,10); 
    $name.=" " x (10-length($name)); 
    $phylipText.="$name  "; 
    if($seenTruncName{$name}++){
      die "ERROR: genome names truncated to 10 characters are not unique! Discovered when looking at $genome[$i].";
    }
    for(my $j=0;$j<@genome;$j++){
      $phylipText.=$m{$genome[$i]}{$genome[$j]}."  ";
    }
    $phylipText.= "\n";
  }
  $phylipText=~s/  $//gm;

  # Make the phylip file.
  my $phylip = "$$settings{tempdir}/distances.phylip"; 
  open(PHYLIP,">",$phylip) or die "ERROR: could not open $phylip for writing: $!";
  my @names=glob("$$settings{tempdir}/msh/*.msh");
  print PHYLIP "    ".scalar(@names)."\n";
  print PHYLIP $phylipText;
  close PHYLIP;

  # Create the tree
  logmsg "Creating tree with fneighbor";
  system("fneighbor -datafile $phylip -outfile $$settings{tempdir}/fneighbor.txt -outtreefile $$settings{tempdir}/fneighbor.dnd -treetype neighbor-joining > /dev/null 2>&1");
  die $! if $?;

  # Read and edit the tree
  open(TREE,"$$settings{tempdir}/fneighbor.dnd") or die $!;
  my @treeLine=<TREE>;
  close TREE;
  chomp(@treeLine);
  my $treeContent=join("",@treeLine);
  $treeContent=~s/^\s+|\s+$//gm; # whitespace and newline trim
  $treeContent.="\n";            # However, give it a newline.

  # Return the edited tree
  return $treeContent;
}

sub usage{
  "$0: use distances from Mash (min-hash algorithm) to make a NJ tree
  Usage: $0 *.fastq.gz > tree.dnd
  --numcpus      1

  MASH SKETCH OPTIONS
  --genomesize   5000000
  --mindepth     2     
  "
}
