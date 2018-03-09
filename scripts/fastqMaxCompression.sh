#!/bin/bash

set -e 

DIR=$1

if [ "$DIR" == "" ]; then
  echo "Finds fastq.gz files and runs max compression on them";
  echo "Usage: $0 dir"
  exit 1
fi

TMP=$(mktemp --directory FASTQMAXCOMPRESSION.XXXXXX --tmpdir=$TMPDIR)
trap "{ rm -rf $TMP; }" EXIT
export TMP

find $DIR -iname '*.fastq.gz' -or -iname '*.fq.gz' | xargs -P 1 -n 1 bash -c '
  if [ "$(file $0 | grep -m 1 -o "max compression" | head -n 1)" != "" ]; then
    echo "Skipping $0 bc it is already max compressed"
    exit 0
  fi

  originalSize=$(du $0 | cut -f 1)
  
  tmpfile=$(mktemp --tmpdir=$TMP MAX.XXXXXX --suffix=.fastq.gz)
  trap "{ rm -f $tmpfile; }" EXIT

  echo "$0 => $tmpfile"
  gzip -dc $0 | gzip -9c > $tmpfile && \
    mv -v $tmpfile $0

  newsize=$(du $0 | cut -f 1);
  savings=$(printf "%0.2f" $(echo "$newsize/$originalSize" | bc -l));

  echo "New file is $savings of original"
'

find $DIR -name '*.fastq' -or '*.fq' | xargs -P 1 -n 1 bash -c '
  tmpfile=$(mktemp --tmpdir=$TMP MAX.XXXXXX --suffix=.fastq.gz)
  trap "{ rm -f $tmpfile; }" EXIT

  gzip -c9 $0 > $tmpfile && \
    mv $tmpfile $0.gz && \
    rm $0
'

