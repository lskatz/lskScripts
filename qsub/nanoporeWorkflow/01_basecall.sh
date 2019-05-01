#!/bin/bash
#$ -o wtdbg2.log
#$ -j y
#$ -N wtdbg2
#$ -pe smp 2-16
#$ -V -cwd
#$ -l gpu=1
set -e

source /etc/profile.d/modules.sh
module purge

NSLOTS=${NSLOTS:=24}

OUTDIR=$1
FAST5DIR=$2
GENOMELENGTH=5000000 # TODO make this a parameter
LONGREADCOVERAGE=50  # How much coverage to target with long reads

set -u

if [ "$FAST5DIR" == "" ]; then
    echo "Usage: $0 outdir fast5dir/"
    echo "  The outdir will represent each sample in a 'barcode__' subdirectory"
    exit 1;
fi;

# Setup any debugging information
date
hostname

# Setup tempdir
tmpdir=$(mktemp -p . -d guppy.XXXXXX)
trap ' { echo "END - $(date)"; rm -rf $tmpdir; } ' EXIT
mkdir $tmpdir/log
echo "$0: temp dir is $tmpdir";

module purge
module load guppy/2.3.5

# Base calling
guppy_basecaller -i $FAST5DIR -s $tmpdir/fastq --gpu_runners_per_device 96 --cpu_threads_per_caller $NSLOTS --num_callers $NSLOTS --flowcell FLO-MIN106 --kit SQK-LSK109 --qscore_filtering 7 --enable_trimming yes --hp_correct yes -r

# Demultiplex.  -r for recursive fastq search.
guppy_barcoder -t $NSLOTS -r -i $tmpdir/fastq -s $tmpdir/demux
# retain the sequencing summary
ln -v $tmpdir/fastq/sequencing_summary.txt $tmpdir/demux

# Make a relative path symlink to the sequencing summary
# file for each barcode subdirectory
for barcodeDir in $tmpdir/demux/barcode[0-9]* $tmpdir/demux/unclassified; do
  ln -sv ../sequencing_summary.txt $barcodeDir/;
done

cp -rv $tmpdir/demux $OUTDIR

