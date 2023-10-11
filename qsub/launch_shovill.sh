#!/bin/bash
# Runs a shovill assembly.  Run with no options for usage.
# Author: Lee Katz <lkatz@cdc.gov>

#$ -S /bin/bash
#$ -pe smp 4-16
#$ -cwd -V
#$ -o shovill.log
#$ -j y
#$ -N Shovill
#$ -q all.q

R1=$1
R2=$2
outdir=$3
NSLOTS=${NSLOTS:=12}

source /etc/profile.d/modules.sh
scriptname=$(basename $0);

if [ "$outdir" == "" ]; then
  echo "Usage: $scriptname reads.1.fastq.gz reads.2.fastq.gz outdir"
  exit 0;
fi;

set -e
set -u

module load shovill

shovill --check

shovill --outdir $outdir --R1 $R1 --R2 $R2 --assembler skesa

