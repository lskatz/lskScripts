#!/bin/bash

SRR=$1

script=$(basename $0);
if [ "$SRR" == "" ]; then
  echo "Downloads a fastq properly using fastq-dump"
  echo "Usage: $script SRR_accession"
  exit 1;
fi

fastq-dump --accession $SRR --outdir . --defline-seq '@$ac.$si/$ri' --defline-qual '+' --split-files --skip-technical --dumpbase --clip

exit $? # carry over the exit code
