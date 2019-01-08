#!/bin/bash
#$ -o wtdbg2.log
#$ -j y
#$ -N wtdbg2
#$ -pe smp 2-16
#$ -V -cwd
set -e

source /etc/profile.d/modules.sh
module purge
module load samtools/1.4.1
module load minimap2

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
    echo "Usage: $0 out.fasta reads.fastq.gz [reads2.fastq.gz...]"
    exit 1;
fi;

which minimap2 wtdbg2 samtools wtpoa-cns

tmpdir=$(mktemp -p . -d wtdbg2.XXXXXX)
trap ' { echo "END - $(date)"; rm -rf $tmpdir; } ' EXIT
mkdir $tmpdir/log

# Combine reads.
# Use zcat -f -- so that it doesn't matter if it's compressed or not
zcat -f -- $READS > $tmpdir/reads.fastq

# Find the desired read length by making a table of
# sorted read lengths vs cumulative coverage.
# First command: find read lengths. Slow step.
LENGTHS=$tmpdir/readlengths.txt
cat $READS | perl -lane 'next if($. % 4 != 2); print length($_);' > $LENGTHS

# Second command: get the table but stop when it gets to
# the desired coverage.
# This is relatively fast.
MINLENGTH=$(sort -rn $LENGTHS | perl -lane 'chomp; $minlength=$_; $cum+=$minlength; $cov=$cum/'$GENOMELENGTH'; last if($cov > '$LONGREADCOVERAGE'); END{print $minlength;}')

echo "Min length for $LONGREADCOVERAGE coverage will be $MINLENGTH";

# Assemble.
wtdbg2 -t $NSLOTS -i $tmpdir/reads.fastq -fo $tmpdir/$PREFIX.wtdbg2 -p 19 -AS 2 -s 0.05 -L $MINLENGTH -g $GENOMELENGTH -X $LONGREADCOVERAGE
# Generate the actual assembly using wtpoa-cns
wtpoa-cns -t $NSLOTS -i $tmpdir/$PREFIX.wtdbg2.ctg.lay.gz -o $tmpdir/$(basename $OUT)

# Polishing the assembly with nanopore reads
# Map reads vs assembly and convert to bam file in one go
minimap2 -t $NSLOTS -x map-ont -a $tmpdir/$(basename $OUT) $tmpdir/reads.fastq | \
  samtools view -Sb - >$tmpdir/prefix.ctg.lay.map.bam 

# polish the assembly with mapped reads using wtpoa-cns
samtools view $tmpdir/prefix.ctg.lay.map.bam | wtpoa-cns -t $NSLOTS -d $tmpdir/$(basename $OUT) -i - -fo $OUT ||\
  cp -v $tmpdir/$(basename $OUT) $OUT


