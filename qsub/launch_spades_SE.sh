#!/bin/bash -l
#$ -S /bin/bash
#$ -pe smp 4-16
#$ -cwd -V
#$ -o spades.log -j y
#$ -N SPAdes3.1.0

reads=$1
out=$2
fasta=$3

if [ "$out" == "" ]; then
  echo "Usage: $0 reads.fastq.gz output/ out.fasta"
  exit 1;
fi;

module load SPAdes/3.1.0
if [ $? -gt 0 ]; then echo "unable to load spades 3.1.0"; exit 1; fi;

NSLOTS=${NSLOTS:=1}

spades.py -s $reads -o $out -t $NSLOTS
if [ $? -gt 0 ]; then echo "problem with spades 3.1.0"; exit 1; fi;

if [ "$fastaOut" != "" ]; then
  cp -v "$out/scaffolds.fasta" $fastaOut
  if [ $? -gt 0 ]; then echo "problem with copying $out/scaffolds.fasta => $fastaOut"; exit 1; fi;
  rm -rf "$out";
  if [ $? -gt 0 ]; then echo "problem with removing the directory $out"; exit 1; fi;
fi

