#!/bin/bash -l

#$ -S /bin/bash
#$ -pe smp 10-16
#$ -cwd -V
#$ -o annotation.log
#$ -j y
#$ -N cgpAnnotate
#$ -q all.q

project=$1

if [ "$project" == "" ]; then
  echo "Usage: " $(basename $0) " CGP-project/"
  exit 1;
fi;

NSLOTS=${NSLOTS:=1}

run_pipeline annotate -p "$project" --numcpus $NSLOTS --skip INTERPRO
if [ $? -gt 0 ]; then echo "ERROR with run_pipeline annotatE"; exit 1; fi;
