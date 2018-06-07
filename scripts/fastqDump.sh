#!/bin/bash

SRR=$1

script=$(basename $0);
if [ "$SRR" == "" ]; then
  echo "Downloads a fastq properly using fastq-dump"
  echo "Usage: $script SRR_accession"
  exit 1;
fi

tempdir=$(mktemp --directory --tmpdir=$TMPDIR $(basename $0).XXXXXX)
trap "{ rm -rf $tempdir; }" EXIT SIGINT SIGTERM
echo "Files will temporarily be stored in $tempdir"

fastq-dump --gzip --accession $SRR --outdir $tempdir --defline-seq '@$ac.$si/$ri' --defline-qual '+' --split-files --skip-technical --dumpbase --clip

if [ $? -gt 0 ]; then 
  echo "ERROR with fastq-dump and $SRR"
  exit 1
fi

mv $tempdir/* .
