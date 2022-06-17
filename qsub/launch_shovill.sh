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

export PATH=$PATH:$HOME/bin/shovill-1.0.9/bin
module load SPAdes/3.13.0 Skesa/2.3.0 megahit/1.1.2 velvet/1.2.10 lighter/1.1.1 flash/1.2.11 samtools/1.9 bwa/0.7.17 Mash/2.0 seqtk/1.3 pilon/1.22 trimmomatic/0.35 perl/5.16.1-MT

shovill --check

#shovill --outdir $outdir --R1 $R1 --R2 $R2
shovill --outdir $outdir --R1 $R1 --R2 $R2 --assembler skesa

