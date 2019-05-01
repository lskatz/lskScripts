#!/bin/bash
#$ -o wtdbg2.log
#$ -j y
#$ -N wtdbg2
#$ -pe smp 2-16
#$ -V -cwd
set -e

source /etc/profile.d/modules.sh
module purge

NSLOTS=${NSLOTS:=24}

FASTQDIR=$1

GENOMELENGTH=5000000 # TODO make this a parameter
LONGREADCOVERAGE=50  # How much coverage to target with long reads

set -u

if [ "$FASTQDIR" == "" ]; then
    echo "Usage: $0 projectdir"
    exit 1;
fi;

# Setup any debugging information
date
hostname

# Setup tempdir
tmpdir=$(mktemp -p . -d prepfastq.XXXXXX)
#trap ' { echo "END - $(date)"; rm -rf $tmpdir; } ' EXIT
mkdir $tmpdir/log
echo "$0: temp dir is $tmpdir";

dir=$FASTQDIR

BARCODE=$(basename $dir)

# Gzip them all
uncompressed=$(\ls $dir/*.fastq 2>/dev/null || true)
if [ "$uncompressed" != "" ]; then
  echo "$uncompressed" | xargs -P $NSLOTS gzip 
fi

# Put all the individual gzip fastqs into a subdir,
# Concatenate them, and then remove them.
# Keep the aggregate fastq file.
mkdir -p $dir/fastqChunks
mv $dir/*.fastq.gz $dir/fastqChunks
cat $dir/fastqChunks/*.fastq.gz > $dir/all.fastq.gz
rm -rf $dir/fastqChunks

# Combine reads and count lengths in one stream
LENGTHS=$dir/readlengths.txt.gz
zcat $dir/all.fastq.gz | perl -lne '
  next if($. % 4 != 2);
  print length($_);
' | sort -rn | gzip -cf > $LENGTHS;

