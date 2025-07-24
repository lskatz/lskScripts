#!/bin/bash

SRR=$1

R1="${SRR}_1.fastq.gz";
R2="${SRR}_2.fastq.gz";
R1uncompressed="${SRR}_1.fastq"
R2uncompressed="${SRR}_2.fastq"

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

# Check if fasten_sort is in the path and if not, quit
echo "Checking dependency paths"
which fasten_sort
which fasterq-dump || which fastq-dump
which gzip
#which esearch efetch

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
  set +x
else
  cd $tempdir
  #fasterq-dump $SRR --print-read-nr --threads 1 --outdir $tempdir --split-files --skip-technical 
  fasterq-dump $SRR --threads 1 --outdir $tempdir --split-files --skip-technical 
  if [ $? -gt 0 ]; then 
    echo "ERROR with fasterq-dump and $SRR"
    exit 1
  fi
  if [ ! -e "$R1uncompressed" ]; then 
    echo "ERROR: R1uncompressed not present in filename $R1uncompressed";
    ls -lhA $tempdir
    exit 1;
  fi
  set +x
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
mv -v $tempdir/$R1uncompressed $tempdir/unsorted_1.fastq
mv -v $tempdir/$R2uncompressed $tempdir/unsorted_2.fastq
cat $tempdir/unsorted_1.fastq $tempdir/unsorted_2.fastq | \
  fasten_shuffle | \
  fasten_sort --sort-by SEQ --paired-end | \
  fasten_progress --print --id "sort-reads_$tempdir" --update-every 100000 | \
  fasten_shuffle -d -1 $tempdir/$R1uncompressed -2 $tempdir/$R2uncompressed

for i in $tempdir/$R1uncompressed $tempdir/$R2uncompressed; do
  if [ -e $i ]; then
    gzip -3 -v $i
  fi
done
rm -vf $tempdir/unsorted_1.fastq $tempdir/unsorted_2.fastq

ls -lhd $tempdir
ls -lh $tempdir/*
for i in $tempdir/{$R1,$R2}; do
  if [ -e $i ]; then
    mv -v $i .
  fi
done


