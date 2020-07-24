#!/bin/bash

module purge
module load sratoolkit/2.9.1

SRR=$1

script=$(basename $0);
if [ "$SRR" == "" ]; then
  echo "Downloads a fastq properly using fastq-dump"
  echo "Usage: $script SRR_accession"
  exit 1;
fi

if [ -e "${SRR}.fastq.gz" ]; then
  echo "${SRR}.fastq.gz is already present."
  exit 1
fi
if [ -e "${SRR}_1.fastq.gz" ]; then
  echo "${SRR}_1.fastq.gz is already present."
  exit 1
fi

tempdir=$(mktemp --directory --tmpdir=$TMPDIR $(basename $0).XXXXXX)
trap "{ rm -rf $tempdir; }" EXIT SIGINT SIGTERM
echo "Files will temporarily be stored in $tempdir"

# Decide whether to run fastq-dump or fasterq-dump
if [ "$(which fasterq-dump 2>/dev/null)" == "" ]; then
  fastq-dump --gzip --accession $SRR --outdir $tempdir --defline-seq '@$ac.$si/$ri' --defline-qual '+' --split-files --skip-technical --dumpbase --clip
  if [ $? -gt 0 ]; then 
    echo "ERROR with fastq-dump and $SRR"
    exit 1
  fi
else
  cd $tempdir
  fasterq-dump $SRR --print-read-nr --threads 1 --outdir $tempdir --split-files --skip-technical 
  if [ $? -gt 0 ]; then 
    echo "ERROR with fasterq-dump and $SRR"
    exit 1
  fi
  cd -

  # Compress fastq files
  for fastq in $tempdir/*.fastq; do 
    # Remove quality defline with perl before compressing
    perl -lane '
      if($. % 4 == 3){
        $_="+";
      }
      print;
    ' < $fastq > $fastq.tmp;
    mv $fastq.tmp $fastq
    gzip $fastq &
  done;
  wait;
fi


mv -v $tempdir/* .
