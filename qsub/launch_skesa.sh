#!/bin/bash -l
#$ -S /bin/bash
#$ -pe smp 1-16
#$ -cwd -V
#$ -o skesa.log -j y
#$ -N Skesa_2.0_2
#$ -q all.q

reads=$1
fastaOut=$2

if [ "$fastaOut" == "" ]; then
  echo "Usage: $0 shuffled.fastq.gz out.fasta"
  exit 1;
fi;

module load Skesa/2.0_2
if [ $? -gt 0 ]; then echo "unable to load Skesa v2"; exit 1; fi;

if [ -e "$fastaOut" ]; then
  echo "$fastaOut already exists. Exiting.";
  exit 1;
fi

NSLOTS=${NSLOTS:=1}
skesa --cores $NSLOTS --fastq $reads --gz --use_paired_ends > $fastaOut

