#!/bin/bash
#$ -o wtdbg2.log
#$ -j y
#$ -N wtdbg2
#$ -pe smp 2-16
#$ -V -cwd
set -e

source /etc/profile.d/modules.sh
module purge
module load tabix samtools/1.3.1 minimap2 # medaka

NSLOTS=${NSLOTS:=1}
#NSLOTS=24

OUT=$1
shift || true
READS=$@
GENOMELENGTH=5000000 # TODO make this a parameter
LONGREADCOVERAGE=50  # How much coverage to target with long reads

set -u

PREFIX=$(basename $OUT .fasta)

if [ "$READS" == "" ]; then
    echo "Usage: $0 out.fasta reads.fastq.gz [reads2.fastq.gz...] [reads.fast5...]"
    exit 1;
fi;

date
hostname
which wtdbg2 wtpoa-cns

tmpdir=$(mktemp -p . -d wtdbg2.XXXXXX)
trap ' { echo "END - $(date)"; rm -rf $tmpdir; } ' EXIT
mkdir $tmpdir/log
mkdir $tmpdir/fast5

# Combine reads.
# Use zcat -f -- so that it doesn't matter if it's compressed or not
# Use fast compression because we value speed here.  Anyways, this 
# is the temp dir and will be cleaned up.
# Any compression in the first place will show some speed
# up with disk reading later on.
for i in $READS; do
  if [[ "$i" =~ fast5$ ]]; then
    ln -v $i $tmpdir/fast5/;
  else
    zcat -v -f -- $READS
  fi
done |\
  gzip -1c > $tmpdir/reads.fastq.gz

# Find the desired read length by making a table of
# sorted read lengths vs cumulative coverage.
# First command: find read lengths. Slow step.
LENGTHS=$tmpdir/readlengths.txt
zcat $tmpdir/reads.fastq.gz | perl -lne 'next if($. % 4 != 2); print length($_);' > $LENGTHS

# Second command: get the table but stop when it gets to
# the desired coverage.
# This is relatively fast.
MINLENGTH=$(sort -rn $LENGTHS | perl -lane 'chomp; $minlength=$_; $cum+=$minlength; $cov=$cum/'$GENOMELENGTH'; last if($cov > '$LONGREADCOVERAGE'); END{print $minlength;}')

echo "Min length for $LONGREADCOVERAGE coverage will be $MINLENGTH";

# Assemble.
wtdbg2 -t $NSLOTS -i $tmpdir/reads.fastq.gz -fo $tmpdir/$PREFIX.wtdbg2 -p 19 -AS 2 -s 0.05 -L $MINLENGTH -g $GENOMELENGTH -X $LONGREADCOVERAGE
# Generate the actual assembly using wtpoa-cns
wtpoa-cns -t $NSLOTS -i $tmpdir/$PREFIX.wtdbg2.ctg.lay.gz -o $tmpdir/$(basename $OUT)

cp -v $tmpdir/$(basename $OUT) $OUT

# Polish
#medaka_consensus -i $tmpdir/reads.fastq.gz -d ${DRAFT} -o ${CONSENSUS} -t ${NPROC}

# Nanopolish is the minion polisher with minion data. I
# combined a few commands here and made a subshell in 
# parentheses.
module purge
module load nanopolish
module load minimap2
module load samtools/1.8
# Index the reads
nanopolish index -d $tmpdir/fast5 $tmpdir/reads.fastq.gz
# Map the reads to get a bam
minimap2 -x map-ont -t $NSLOTS $tmpdir/$(basename $OUT) $tmpdir/reads.fastq.gz | \
  samtools view -bS -T $tmpdir/$(basename $OUT) - - > $tmpdir/unsorted.bam
samtools sort -l 1 $tmpdir/unsorted.bam > $tmpdir/reads.bam

# Start a loop based on suggested ranges using nanopolish_makerange.py
# but invoke it with python
for window in $(python $(which nanopolish_makerange.py) $tmpdir/$(basename $OUT)); do
  echo "WINDOW: $window"; 
  nanopolish variants --consensus $tmpdir/consensus.$window.fasta -r $tmpdir/reads.fastq.gz -b $tmpdir/reads.bam -g $tmpdir/$(basename $OUT) -t $NSLOTS --min-candidate-frequency 0.1 --min-candidate-depth 20 -w "$window" > $tmpdir/consensus.$window.vcf; 
done;
nanopolish_merge.py $tmpdir/consensus.*.fasta > $tmpdir/out.polished.fasta
mv -v $tmpdir/out.polished.fasta $OUT

