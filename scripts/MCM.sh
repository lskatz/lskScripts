#!/bin/sh

# runs mauve contig mover
# Author: Lee Katz

ref=$1
draft=$2
output=$3
jar=`dirname $0`/Mauve.jar

if [ "$output" = "" ]; then
  echo "Usage: $0 ref.fasta draft.fasta outputDir";
  echo "NOTE: this script must be in the same directory as Mauve.jar"
  exit 1;
fi

java -Xmx500m -cp $jar org.gel.mauve.contigs.ContigOrderer -output $output -ref $ref -draft $draft
