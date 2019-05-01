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

INDIR=$1
GENOMELENGTH=5000000 # TODO make this a parameter
LONGREADCOVERAGE=50  # How much coverage to target with long reads

set -u

if [ "$INDIR" == "" ]; then
    echo "Usage: $0 projectDir"
    exit 1;
fi;

# Setup any debugging information
date
hostname

# Setup tempdir
tmpdir=$(mktemp -p . -d wtdbg2.XXXXXX)
trap ' { echo "END - $(date)"; rm -rf $tmpdir; } ' EXIT
mkdir $tmpdir/log
echo "$0: temp dir is $tmpdir";

dir=$INDIR

LENGTHS=$dir/readlengths.txt
# Determine minimum read length for desired coverage.
# Do this by reading lengths from biggest to smallest,
# stopping when we get to the desired coverage and saving
# that read length.
MINLENGTH=$(zcat "$LENGTHS" | sort -rn | perl -lane 'chomp; $minlength=$_; $cum+=$minlength; $cov=$cum/'$GENOMELENGTH'; last if($cov > '$LONGREADCOVERAGE'); END{print $minlength;}')
echo "Min length for $LONGREADCOVERAGE coverage will be $MINLENGTH";

# Assemble.
module purge
module load wtdbg2/2.4
wtdbg2 -t $NSLOTS -i $dir/all.fastq.gz -fo $dir/wtdbg2 -p 19 -AS 2 -s 0.05 -L $MINLENGTH -g $GENOMELENGTH -X $LONGREADCOVERAGE

# Generate the actual assembly using wtpoa-cns
wtpoa-cns -t $NSLOTS -i $dir/wtdbg2.ctg.lay.gz -o $dir/unpolished.fasta 

