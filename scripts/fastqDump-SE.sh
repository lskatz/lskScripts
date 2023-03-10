#!/bin/bash

SRR=$1

R1="${SRR}.fastq.gz";
R1uncompressed="${SRR}.fastq"

script=$(basename $0);
if [ "$SRR" == "" ]; then
  echo "Downloads a fastq properly using fastq-dump"
  echo "  Sorts the reads for maximum compression using fasten_sort."
  echo "Usage: $script SRR_accession"
  exit 1;
fi

set -e
set -u

if [ -e "${SRR}.fastq.gz" ]; then
  echo "${SRR}.fastq.gz is already present."
  exit 1
fi
if [ -e "$R1" ]; then
  echo "$R1 is already present."
  exit 1
fi

module purge
#module load sratoolkit/2.9.1
module load sratoolkit/2.11.3

# Check if fasten_sort is in the path and if not, quit
echo "Checking dependency paths"
which fasten_sort
which fasterq-dump || which fastq-dump

tempdir=$(mktemp --directory --tmpdir=$TMPDIR $(basename $0).XXXXXX)
trap "{ rm -rf $tempdir; }" EXIT SIGINT SIGTERM
echo "Files will temporarily be stored in $tempdir"

# Decide whether to run fastq-dump or fasterq-dump
fasterqDump="$(which fasterq-dump 2>/dev/null)";
if [ "$fasterqDump" == "" ]; then
  fastq-dump --accession $SRR --outdir $tempdir --defline-seq '@$ac.$si/$ri' --defline-qual '+' --split-files --skip-technical --dumpbase --clip
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
  done;
fi

# Intense compression
mv -v $tempdir/$R1uncompressed $tempdir/unsorted.fastq
cat $tempdir/unsorted.fastq | \
  fasten_sort --sort-by SEQ | \
  fasten_progress --print --id sort-reads --update-every 100000 | \
  gzip -vc9 > $tempdir/$R1uncompressed.gz

rm -v $tempdir/unsorted.fastq

ls -lhd $tempdir
ls -lh $tempdir/*
mv -v $tempdir/$R1 .

