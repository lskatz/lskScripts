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
  GetOptions($settings,qw(help numcpus=i genomesize=i mindepth=i truncLength=i warn-on-duplicate)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{genomesize}||=5000000;
  $$settings{mindepth}||=2;
  $$settings{truncLength}||=10;  # how long a genome name is
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

  logmsg "$0 on ".scalar(@reads)." files";

  validateFastq(\@reads,$settings);

  my $mshList=sketchAll(\@reads,$settings);

  my $distances=mashDistance($mshList,$settings);

  my $treeContent = distancesToTree($distances,$settings);

  print $treeContent;

  return 0;
}

sub validateFastq{
  my($reads,$settings)=@_;
  
  my %seen;
  for my $r(@$reads){
    my $trunc=_truncateFilename($r,$settings);
    if($seen{$trunc}){
      my $msg="I have already seen $r as $seen{$trunc} (truncated name: $trunc)";
      if($$settings{'warn-on-duplicate'}){
        logmsg "WARNING: $msg";
      } else {
        die "ERROR: $msg";
      }
    }
    $seen{$trunc}=$r;
  }

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

  my $mashlist="$$settings{tempdir}/msh.txt";
  open(MASHLIST,">",$mashlist) or die "ERROR: could not open $mashlist for writing: $!";
  for(@thr){
    my $mashfiles=$_->join;
    for my $file(@$mashfiles){
      print MASHLIST $file."\n";
    }
  }
  close MASHLIST;

  return $mashlist;
}

# Individual mash sketch
sub mashSketch{
  my($sketchDir,$Q,$settings)=@_;

  my @msh;
  while(defined(my $fastq=$Q->dequeue)){
    logmsg "Sketching $fastq";
    my $outPrefix="$sketchDir/".basename($fastq);
    if(-e "$outPrefix.msh"){
      logmsg "WARNING: ".basename($fastq)." was already mashed. You need unique filenames for this script. This file will be skipped: $fastq";
      next;
    }
    if(-s $fastq < 1){
      logmsg "WARNING: $fastq is a zero byte file. Skipping.";
      next;
    }
    system("mash sketch -k 21 -s 10000 -m $$settings{mindepth} -c 10 -g $$settings{genomesize} -o $outPrefix $fastq > /dev/null 2>&1");
    die if $?;

    push(@msh,"$outPrefix.msh");
  }

  return \@msh;
}

# Parallelized mash distance
sub mashDistance{
  my($mshList,$settings)=@_;

  open(MASHLIST,"<",$mshList) or die "ERROR: could not open $mshList for reading: $!";
  my @msh=<MASHLIST>; chomp(@msh);
  close MASHLIST;

  my $outdir="$$settings{tempdir}/dist";
  mkdir($outdir);

  my $mshQueue=Thread::Queue->new(@msh);
  my @thr;
  for(0..$$settings{numcpus}-1){
    $thr[$_]=threads->new(\&mashDist,$outdir,$mshQueue,$mshList,$settings);
  }

  $mshQueue->enqueue(undef) for(@thr);

  my $distfile="$$settings{tempdir}/distances.tsv";
  open(DIST,">",$distfile) or die "ERROR: could not open $distfile for writing: $!";
  for(@thr){
    my $distfiles=$_->join;
    for my $file(@$distfiles){
      # Print the contents of each dist file to the
      # main dist file.
      open(ONEDISTFILE,"<",$file) or die "ERROR: could not open $file for reading: $!";
      while(<ONEDISTFILE>){
        print DIST $_;
      }
      close ONEDISTFILE;
    }
  }
  close DIST;

  return $distfile;
}

# Individual mash distance
sub mashDist{
  my($outdir,$mshQueue,$mshList,$settings)=@_;
  my @dist;
  while(defined(my $msh=$mshQueue->dequeue)){
    my $outfile="$outdir/".basename($msh).".tsv";
    logmsg "Distances for $msh";
    system("mash dist -t $msh -l $mshList > $outfile");
    die if $?;

    push(@dist,$outfile);
  }

  return \@dist;
}

# 1. Read the mash distances
# 2. Create a phylip file
# 3. Run emboss's fneighbor
sub distancesToTree{
  my($distances,$settings)=@_;

  logmsg "Reading the distances file at $distances";
  open(MASHDIST,"<",$distances) or die "ERROR: could not open $distances for reading: $!";

  my $id="UNKNOWN"; # Default ID in case anything goes wrong
  my %m; #matrix for distances
  while(<MASHDIST>){
    chomp;
    if(/^#query\s+(.+)/){
      #$id=basename($1,@fastqExt);
      $id=_truncateFilename($1,$settings);
    } else {
      my @F=split(/\t/,$_);
      #$F[0]=basename($F[0],@fastqExt);
      $F[0]=_truncateFilename($F[0],$settings);
      $m{$id}{$F[0]}=sprintf("%0.6f",$F[1]);
    }
  }
  close MASHDIST;

  # Create the phylip file.
  # Make the text first so that we can edit it a bit.
  logmsg "Creating the distance matrix file for fneighbor.";
  my %seenTruncName;
  my $phylipText="";
  my @genome=sort{$a cmp $b} keys(%m);
  for(my $i=0;$i<@genome;$i++){ 
    my $name=_truncateFilename($genome[$i],$settings);
    $phylipText.="$name  "; 
    if($seenTruncName{$name}++){
      #die "ERROR: genome names truncated to $$settings{truncLength} characters are not unique! Discovered when looking at $genome[$i].";
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
  print PHYLIP "    ".scalar(@genome)."\n";
  print PHYLIP $phylipText;
  close PHYLIP;

  # Create the tree
  logmsg "Creating tree with fneighbor";
  system("fneighbor -datafile $phylip -outfile $$settings{tempdir}/fneighbor.txt -outtreefile $$settings{tempdir}/fneighbor.dnd -treetype neighbor-joining 1>&2");
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

#######
# Utils
#######

sub _truncateFilename{
  my($file,$settings)=@_;
  my $name=basename($file,@fastqExt);
  $name=substr($name,0,$$settings{truncLength}); 
  $name.=" " x ($$settings{truncLength}-length($name)); 
  return $name;
}

sub usage{
  "$0: use distances from Mash (min-hash algorithm) to make a NJ tree
  Usage: $0 *.fastq.gz > tree.dnd
  --numcpus            1
  --truncLength        10   How many characters to keep from a filename
  --warn-on-duplicate       Warn instead of die when a duplicate
                            genome name is found

  MASH SKETCH OPTIONS
  --genomesize   5000000
  --mindepth     2     
  "
}
