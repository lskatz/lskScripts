#!/bin/bash -l
#$ -S /bin/bash
#$ -pe smp 4-16
#$ -cwd -V
#$ -o spades.log
#$ -j y
#$ -N SPAdes3.9.0
#$ -q all.q

fwd=$1
rev=$2
out=$3
fastaOut=$4

if [ "$out" == "" ]; then
  echo "Usage: $0 fwd.fastq.gz rev.fastq.gz output/ [out.fasta]"
  echo "  if out.fasta is given, then the output directory will be removed and the scaffolds.fasta file will be saved"
  exit 1;
fi;

module load SPAdes/3.9.0
if [ $? -gt 0 ]; then echo "unable to load spades 3.9.0"; exit 1; fi;
module load Python/2.7.12

NSLOTS=${NSLOTS:=1}
spades.py -1 $fwd -2 $rev --careful -o $out -t $NSLOTS 
if [ $? -gt 0 ]; then echo "problem with spades 3.9.0"; exit 1; fi;

if [ "$fastaOut" != "" ]; then
  cp -v "$out/scaffolds.fasta" $fastaOut
  if [ $? -gt 0 ]; then echo "problem with copying $out/scaffolds.fasta => $fastaOut"; exit 1; fi;
  rm -rf "$out";
  if [ $? -gt 0 ]; then echo "problem with removing the directory $out"; exit 1; fi;
fi
