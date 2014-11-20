#!/bin/bash -l
#$ -pe smp 16
#$ -cwd -V
#$ -o velvet.log -j y
#$ -N velvet

module load velvet/1.2.10;
if [ $? -gt 0 ]; then echo "unable to load velvet/1.2.10"; exit 1; fi;

reads=$1
out=$2

# number of cpus is either set by SGE, or is just 1
NSLOTS=${NSLOTS:=1}
echo $NSLOTS

if [ "$out" == "" ]; then
  echo "Usage: $0 reads.fastq.gz output/"
  exit 1;
fi;

VelvetOptimiser.pl -s 55 -e 99 -d $out -p $out -t $NSLOTS -f "-fastq.gz -shortPaired $reads"
if [ $? -gt 0 ]; then echo "problem with VelvetOptimiser"; exit 1; fi;
