#!/bin/bash -l

#$ -S /bin/bash
#$ -pe smp 10-16
#$ -cwd -V
#$ -o trimClean.log -j y
#$ -N cgpTrimClean
#$ -q all.q

reads=$1
out=$2

if [ "$out" == "" ]; then
  echo "Usage: $0 reads.fastq.gz cleaned.fastq.gz"
  exit 1;
fi;

NSLOTS=${NSLOTS:=1}

run_assembly_trimClean.pl -i $reads -o $out --auto --nosingletons --numcpus $NSLOTS
if [ $? -gt 0 ]; then echo "ERROR with run_assembly_trimClean.pl"; exit 1; fi;
