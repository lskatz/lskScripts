#!/bin/bash
#$ -o wtdbg2.log
#$ -j y
#$ -N wtdbg2
#$ -pe smp 2-16
#$ -V -cwd
set -e

source /etc/profile.d/modules.sh
module purge

NSLOTS=${NSLOTS:=1}
NSLOTS=36

OUTDIR=$1
shift || true
FAST5DIR=$1
GENOMELENGTH=5000000 # TODO make this a parameter
LONGREADCOVERAGE=50  # How much coverage to target with long reads

set -u

if [ "$FAST5DIR" == "" ]; then
    echo "Usage: $0 fastaprefix fast5dir/"
    exit 1;
fi;

# Setup any debugging information
date
hostname

# Setup tempdir
tmpdir=$(mktemp -p . -d wtdbg2.XXXXXX)
#trap ' { echo "END - $(date)"; rm -rf $tmpdir; } ' EXIT
mkdir $tmpdir/log
echo "$0: temp dir is $tmpdir";

module purge
module load guppy/2.3.5

# Base calling
guppy_basecaller -i $FAST5DIR -s $tmpdir/fastq --cpu_threads_per_caller 1 --num_callers $NSLOTS --flowcell FLO-MIN106 --kit SQK-LSK109 --qscore_filtering 7 --enable_trimming yes --hp_correct yes -r

# Demultiplex.  -r for recursive fastq search.
guppy_barcoder -t $NSLOTS -r -i $tmpdir/fastq -s $tmpdir/demux

# Each demultiplexed set of reads is in a subfolder called barcode__
for dir in $tmpdir/demux/barcode*; do
  if [ ! -d $dir ]; then continue; fi;

  BARCODE=$(basename $dir)

  echo "$0: Let's look at $dir";
  gzip -1c $dir/*.fastq > $dir/all.fastq.gz && \
    rm $dir/*.fastq

  # Combine reads and count lengths in one stream
  LENGTHS=$dir/readlengths.txt
  zcat $dir/all.fastq.gz | perl -lne '
    next if($. % 4 != 2);
    print length($_);
  ' > $LENGTHS;

  # Determine minimum read length for desired coverage.
  # Do this by reading lengths from biggest to smallest,
  # stopping when we get to the desired coverage and saving
  # that read length.
  MINLENGTH=$(sort -rn $LENGTHS | perl -lane 'chomp; $minlength=$_; $cum+=$minlength; $cov=$cum/'$GENOMELENGTH'; last if($cov > '$LONGREADCOVERAGE'); END{print $minlength;}')
  echo "Min length for $LONGREADCOVERAGE coverage will be $MINLENGTH";

  # Assemble.
  module purge
  module load wtdbg2/2.4
  wtdbg2 -t $NSLOTS -i $dir/all.fastq.gz -fo $dir/wtdbg2 -p 19 -AS 2 -s 0.05 -L $MINLENGTH -g $GENOMELENGTH -X $LONGREADCOVERAGE || \
    continue
  # Generate the actual assembly using wtpoa-cns
  wtpoa-cns -t $NSLOTS -i $dir/wtdbg2.ctg.lay.gz -o $dir/unpolished.fasta ||\
    continue

  module purge
  module load nanopolish/0.11.1
  module load minimap2
  module load samtools/1.8
  module load tabix/0.2.6
  # Index the reads
  nanopolish index -d $FAST5DIR $dir/all.fastq.gz
  # Map the reads to get a bam
  minimap2 -a -x map-ont -t $NSLOTS $dir/unpolished.fasta $dir/all.fastq.gz | \
    samtools view -bS -T $dir/unpolished.fasta > $dir/unsorted.bam
  samtools sort -l 1 $dir/unsorted.bam > $dir/reads.bam
  rm $dir/unsorted.bam

  # Start a loop based on suggested ranges using nanopolish_makerange.py
  # but invoke it with python
  python nanopolish_makerange.py $dir/unpolished.fasta | \
    xargs -P $NSLOTS -n 1 bash -c '
      window="$0";
      dir="'$dir'";
      nanopolish variants --consensus -r $dir/all.fastq.gz -b $dir/reads.bam -g $dir/unpolished.fasta -t 1 --min-candidate-frequency 0.1 --min-candidate-depth 20 -w "$window" > $dir/consensus.$window.vcf;
      #bgzip $dir/consensus.$window.vcf # help with some level of compression in the folder
    '
  # Run merge on vcf files
  nanopolish vcf2fasta -g $dir/unpolished.fasta $dir/consensus.*.vcf > $dir/polished.fasta

  cp -v $dir/polished.fasta $OUTDIR/$BARCODE.fasta

done;

