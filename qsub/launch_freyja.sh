#!/bin/bash
# Runs Frejya wastewater pipeline for SARS-CoV-2
# Author: Lee Katz <lkatz@cdc.gov>

#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o frejya
#$ -j y
#$ -N frejya

ref=$1; shift
outdir=$1; shift
reads=$@
NSLOTS=${NSLOTS:=1}

source /etc/profile.d/modules.sh
scriptname=$(basename $0);

if [ "$reads" == "" ]; then
  echo "Usage: $scriptname ref.fasta outdir *_1.fastq[.gz]"
  echo "  R2 reads will be detected automatically based on matchiing filenames"
  exit 0;
fi;

set -e
set -u

module purge

which ivar
which freyja
which samtools
which bowtie2

tmpdir=$(mktemp --directory --suffix=$(basename $0));
trap ' { rm -rf $tmpdir; } ' EXIT

mkdir -v $outdir || echo "WARNING: outdir already exists: $outdir"

for R1 in $reads; do
  name=$(basename $R1 .fastq.gz | perl -plane 's/_\d|\.fastq|\.gz//g');
  sampledir="$tmpdir/$name.freyja"
  echo "START $name in $sampledir"
  mkdir -v $sampledir

  # Get the file extension, if it's .gz
  ext=${R1##*.}

  # Get fastq files into our local tmp folder
  R2=${R1/_1.f/_2.f}

  echo "$R1 $R2";

  cp -v $R1 $R2 $sampledir/
  R1=$sampledir/$(basename $R1)
  R2=$sampledir/$(basename $R2)

  bowtie2 -x $ref -1 $R1 -2 $R2 | samtools view -bhS | samtools sort > $sampledir/sorted.bam

  # Trim primers
  ivar trim -i $sampledir/sorted.bam -p $sampledir/ivar.unsorted
  samtools sort $sampledir/ivar.unsorted.bam > $sampledir/ivar.bam

  # call variants
  samtools mpileup -aa -A -d 600000 -B -Q 0 $sampledir/ivar.bam | ivar variants -p $sampledir/ivar -q 20 -t 0.03 -r $ref

  # abundances
  freyja variants $sampledir/ivar.bam --variants $sampledir/freyja.variants --depths $sampledir/freyja.depths --ref $ref
  freyja demix $sampledir/freyja.variants.tsv $sampledir/freyja.depths --output $sampledir/freyja.demix

  # Clear results we don't need
  (cd $sampledir && rm -v *.bam* *.fastq.gz)
  # Save results
  rsync -av $sampledir $outdir/
done

