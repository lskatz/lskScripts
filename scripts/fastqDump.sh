#!/bin/bash

SRR=$1

R1="${SRR}_1.fastq.gz";
R1uncompressed="${SRR}_1.fastq"

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
if [ -e "$R1" ]; then
  echo "$R1 is already present."
  exit 1
fi

module purge
module load sratoolkit/2.9.1

tempdir=$(mktemp --directory --tmpdir=$TMPDIR $(basename $0).XXXXXX)
trap "{ rm -rf $tempdir; }" EXIT SIGINT SIGTERM
echo "Files will temporarily be stored in $tempdir"

# Decide whether to run fastq-dump or fasterq-dump
if [ "$(which fasterq-dump 2>/dev/null)" == "" ]; then
  set -x
  fastq-dump --gzip --accession $SRR --outdir $tempdir --defline-seq '@$ac.$si/$ri' --defline-qual '+' --split-files --skip-technical --dumpbase --clip
  set +x
  if [ $? -gt 0 ]; then 
    echo "ERROR with fastq-dump and $SRR"
    exit 1
  fi
else
  set -x
  cd $tempdir
  fasterq-dump $SRR --print-read-nr --threads 1 --outdir $tempdir --split-files --skip-technical 
  set +x
  if [ $? -gt 0 ]; then 
    echo "ERROR with fasterq-dump and $SRR"
    exit 1
  fi
  if [ ! -e "$R1uncompressed" ]; then 
    echo "ERROR: R1uncompressed not present in filename $R1uncompressed";
    ls -lhA $tempdir
    exit 1;
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

