#!/bin/bash -l

#$ -S /bin/bash
#$ -pe smp 10-16
#$ -cwd -V
#$ -o prediction.log
#$ -j y
#$ -N cgpPredict
#$ -q all.q

project=$1

if [ "$project" == "" ]; then
  echo "Usage: " $(basename $0) " CGP-project/"
  exit 1;
fi;

NSLOTS=${NSLOTS:=1}

run_pipeline predict -p "$project" --numcpus $NSLOTS
if [ $? -gt 0 ]; then echo "ERROR with run_pipeline predict"; exit 1; fi;
