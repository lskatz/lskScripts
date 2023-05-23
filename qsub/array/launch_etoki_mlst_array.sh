#!/bin/bash -l
# Runs EToKi MLST
# Author: Lee Katz <lkatz@cdc.gov>

#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o EToKi.log
#$ -j y
#$ -N EToKi

outdir=$1
refs=$2
db=$3
shift;shift;
asm=$@

#NSLOTS=${NSLOTS:=1}

source /etc/profile.d/modules.sh
scriptname=$(basename $0);


if [ "$asm" == "" ]; then
  echo "Usage: $scriptname outdir refs.fasta etoki.csv *.fasta"
  exit 0;
fi;

conda activate etoki || echo "could not activate etoki env"

set -e
set -u

which EToKi.py
#EToKi.py configure

tmpdir=$(mktemp --tmpdir=. --directory --suffix=.$(basename $0));
#trap ' { rm -rf $tmpdir; } ' EXIT
mkdir -p $tmpdir/log
mkdir -p $tmpdir/scratch
echo "tmp dir is $tmpdir"

CTRL_FILE="$tmpdir/array.txt"
echo "$asm" | tr ' ' '\n' > $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"

if [ -d "$outdir" ]; then
  echo "WARNING: outdir already exists: $outdir"
  echo "  pausing 2 seconds in case you want to cancel.";
  sleep 2;
fi
mkdir -pv $outdir

#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o EToKi.log
#$ -j y
#$ -N EToKi

qsub -N EToKi -o $tmpdir/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "outdir=$outdir" -v "refs=$refs" -v "db=$db" -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash -l

  set -eu
  hostname

  which EToKi.py
  EToKi.py configure || true
  echo

  asm=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE);
  mkdir /scratch/$USER || true

  samplename=$(basename $asm .fasta)
  b=$(basename $asm)
  echo "Sample name will be $samplename"

  EToKi.py MLSType -i $asm -r $refs -k $samplename -d $db -o $tmpdir/scratch/$b
  mv -v $tmpdir/scratch/$b $outdir/$b

END_OF_SCRIPT

