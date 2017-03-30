#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use Getopt::Long;

use Bio::SeqIO;

$0=basename $0;
sub logmsg { print STDERR "$0: @_\n";}

exit main();

sub main{
  my $settings={};
  GetOptions($settings, qw(query=s ref=s lt=f tempdir=s numcpus=i mauveJar=s help)) or die $!;
  $$settings{numcpus}||=1;
  $$settings{lt}||=0.3;
  $$settings{lt} *= 1000000; # Mb to bp conversion
  $$settings{tempdir}||=tempdir("subtractContigs.XXXXXX", TMPDIR=>1, CLEANUP=>1);

  my $contigsQuery=$$settings{query};
  my $contigsRef=$$settings{ref};
  
  die usage() if(!$contigsQuery || !$contigsRef || $$settings{help});

  my $querySeq=mauve($contigsQuery, $contigsRef, $settings);

  my $seqout=Bio::SeqIO->new(-format=>"fasta");
  $seqout->write_seq(values(%$querySeq));

  return 0;
}

sub mauve{
  my($contigsQuery, $contigsRef, $settings)=@_;
  
  # Detect where the Mauve.jar file is
  $$settings{mauveJar}||=findMauveJar($settings);

  my $tmpdir=tempdir("$$settings{tempdir}/mauveContigMover.XXXXXX");

  my $reference="$$settings{tempdir}/B.fasta";
  my $query    ="$$settings{tempdir}/A.fasta";
  filterContigs($contigsQuery, $query    , $settings);
  filterContigs($contigsRef, $reference, $settings);
  
  system("java -Xmx500m -cp $$settings{mauveJar} org.gel.mauve.contigs.ContigOrderer -output $tmpdir -ref $reference -draft $query --backbone-output='backbone' >&2 ");
  die if $?;

  # Get the last iteration
  my @iterDir=sort {
    my $iterA=$a;
    my $iterB=$b;
    $iterA=~s/\D//g;
    $iterB=~s/\D//g;
    return $iterB <=> $iterA;
  } glob("$tmpdir/*");
  my $lastIter=$iterDir[0];
  my $alnPrefix=basename($lastIter);

  # Read the tab file to find contig starts/stops
  my $tabFile="$lastIter/A_contigs.tab";
  my %coordinates;
  open(my $fh, $tabFile) or die "ERROR: could not read $tabFile: $!";
  while(<$fh>){
    next if(!/^contig/);
    chomp;
    my(undef, $seqid, undef, $strand, $left, $right)=split /\t/;

    # Just save the first instance of the contig
    $coordinates{$seqid}//=[$left,$right];
  }
  close $fh;
  my @querySeqid=sort {$coordinates{$a}[0] <=> $coordinates{$b}[0]} keys(%coordinates);

  my %queryseq;
  my $queryseqin=Bio::SeqIO->new(-file=>"$lastIter/A.fasta");
  while(my $seq=$queryseqin->next_seq){
    $queryseq{$seq->id}=$seq;
  }
  $queryseqin->close;

  #for(@sortedQueryContigs){
  #  print($_."\t".$queryseqLength{$_}."\n");
  #}
  #die;

  # Read the backbone file to find which contigs are
  # aligned or not
  my $backbone="$lastIter/$alnPrefix.backbone";
  my $refseqIndex=0;
  my $queryseqIndex=0;
  open(my $fp, $backbone) or die "ERROR: could not read $backbone: $!";
  my $header=<$fp>; 
  chomp($header);
  my @backboneHeader=split(/\t/, $header);
  while(<$fp>){
    s/^\s+|\s+$//g;
    my($lcbSeq0Left,$lcbSeq0Right,$lcbSeq1Left,$lcbSeq1Right)=split(/\t/);
    
    # Don't filter anything out if there are zeros
    # for the query or ref contigs.
    # Seq0 is the reference; seq1 is the query.
    my $refLength  =abs($lcbSeq0Left-$lcbSeq0Right);
    my $queryLength=abs($lcbSeq1Left-$lcbSeq1Right);
    if($refLength==0 || $queryLength==0){
      next;
    }

    # Which contigs do we exclude based on these coordinates?
    my $queryId="UNKNOWN";
    for(my $i=0;$i<@querySeqid;$i++){
      #print join("\t",$lcbSeq1Left,$lcbSeq1Right,$coordinates{$querySeqid[$i]}[0],$coordinates{$querySeqid[$i]}[1],$querySeqid[$i])."\n";
      if(
        ( # if the LCB is within the contig
          # start is <= than the backbone seq1 left end
          ($coordinates{$querySeqid[$i]}[0] <= $lcbSeq1Left) &&
          # stop is >= than the backbone seq1 right end
          ($coordinates{$querySeqid[$i]}[1] >= $lcbSeq1Right)
        ) ||
        (
          ($coordinates{$querySeqid[$i]}[0] <= $lcbSeq1Left) &&
          ($coordinates{$querySeqid[$i]}[1] >= $lcbSeq1Left)
        ) ||
        (
          ($coordinates{$querySeqid[$i]}[0] <= $lcbSeq1Right) &&
          ($coordinates{$querySeqid[$i]}[1] >= $lcbSeq1Right)
        ) ||
        # If the contig stop coordinate is the same as the LCB stop
        ($coordinates{$querySeqid[$i]}[1] == $lcbSeq1Right) ||
        # If the start is the same as the LCB start
        ($coordinates{$querySeqid[$i]}[0] == $lcbSeq1Left)
      ){
        #print ":)\n";
        logmsg "Delete $querySeqid[$i]";
        delete($queryseq{$querySeqid[$i]});
        #last;
      }
    }
  }
  close $fp;

  return \%queryseq;
}

sub filterContigs{
  my($source,$destination, $settings)=@_;

  my $seqin=Bio::SeqIO->new(-file=>$source);
  my $seqout=Bio::SeqIO->new(-file=>">$destination");
  while(my $seq=$seqin->next_seq){
    next if($seq->length >= $$settings{lt});
    $seqout->write_seq($seq);
  }
}

sub findMauveJar{
  my($settings)=@_;

  $ENV{CLASSPATH}//="";
  for my $path("", split(/:/, $ENV{PATH}.":".$ENV{CLASSPATH})){
    if(-e $path && -f "$path/Mauve.jar"){
      return "$path/Mauve.jar";
    }
  }

  die "ERROR: could not find Mauve.jar in PATH or CLASSPATH, and it was not supplied with --mauveJar";
}

sub usage{
  "$0: subtracts contigs found in b.fasta from a.fasta.

  Usage: $0 -query query.fasta -ref ref.fasta > out.fasta
  All contigs found in ref.fasta will be filtered from query.fasta

  --numcpus    1
  --tempdir    ''     One will be made by default
  --mauveJar   ''     This script will look for Mauve.jar
                      in PATH and CLASSPATH, but you can
                      supply it instead with this arg.
  --lt         0.3    Don't compare contigs greater than
                      this many megabases. Do not output
                      contigs greater than this many
                      megabases.
  "
}
