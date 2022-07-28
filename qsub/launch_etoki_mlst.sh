#!/bin/bash
# Runs EToKi MLST
# Author: Lee Katz <lkatz@cdc.gov>

#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o EToKi.log
#$ -j y
#$ -N EToKi

out=$1
refs=$2
db=$3
asm=$4

#NSLOTS=${NSLOTS:=1}

source /etc/profile.d/modules.sh
scriptname=$(basename $0);


if [ "$asm" == "" ]; then
  echo "Usage: $scriptname out.etoki.fasta refs.fasta etoki.csv assembly.fasta"
  exit 0;
fi;

set -e
set -u

module purge

source ~/.bash_conda > /dev/null 2>&1 || echo "Could not find bash file for loading conda"
conda activate etoki || echo "could not activate etoki env"

EToKi.py configure

tmpdir=$(mktemp --directory --suffix=$(basename $0));
trap ' { rm -rf $tmpdir; } ' EXIT

#mkdir -v $outdir || echo "WARNING: outdir already exists: $outdir"

if [ -e $out ]; then
  echo "ERROR: output file already exists: $out"
  exit 1
fi

samplename=$(basename $asm .fasta)
echo "Sample name will be $samplename"

EToKi.py MLSType -i $asm -r $refs -k $samplename -d $db -o $tmpdir/$(basename $out)
mv -v $tmpdir/$(basename $out) $out
