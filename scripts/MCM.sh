#!/bin/sh

# runs mauve contig mover
# Author: Lee Katz

set -e

ref=$1
draft=$2
output=$3
MAUVE=$(which Mauve)
jar=`dirname $MAUVE`/Mauve.jar
echo "JAR: $jar"

if [ "$output" = "" ]; then
  echo "Usage: $0 ref.fasta draft.fasta outputDir";
  exit 1;
fi

java -Xmx500m -cp $jar org.gel.mauve.contigs.ContigOrderer -output $output -ref $ref -draft $draft

cd $output
lastAlignment=$(\ls -t alignment*/alignment* | grep -P '\d$' | head -n 1)
ln -s $lastAlignment ./alignment.xmfa
cd ..
echo "Best alignment can be found at $output/alignment.xmfa"
