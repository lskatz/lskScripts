#!/bin/bash
#$ -S /bin/bash
#$ -pe smp 12-16
#$ -cwd -V
#$ -o spades.log -j y
#$ -N SPAdes3.1.0

reads=$1
out=$2

if [ "$out" == "" ]; then
  echo "Usage: $0 reads.fastq.gz output/"
  exit 1;
fi;

if [ $? -gt 0 ]; then echo "unable to load spades 3.1.0"; exit 1; fi;

spades.py -s $reads -o $out -t $NSLOTS
if [ $? -gt 0 ]; then echo "problem with spades 3.1.0"; exit 1; fi;
