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

conda env list | grep etoki > /dev/null || \
  echo "DIRE WARNING: 'etoki' virtual env not found"

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
for i in $asm; do
  if [ -e "$outdir/$(basename $i)" ]; then
    continue;
  fi
  echo $i
done > $CTRL_FILE

echo "CTRL_FILE is $CTRL_FILE"

if [ -d "$outdir" ]; then
  echo "WARNING: outdir already exists: $outdir"
  echo "  pausing 2 seconds in case you want to cancel.";
  sleep 2;
fi
mkdir -pv $outdir

qsub -N EToKi -o $tmpdir/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "outdir=$outdir" -v "refs=$refs" -v "db=$db" -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash -l

  set -eu
  hostname

  which EToKi.py
  EToKi.py configure || true
  echo

  tmpdir=$(mktemp --directory $(basename $0).TASK_ID$SGE_TASK_ID.XXXXXX --tmpdir=$TMPDIR)
  trap ' { rm -rf $tmpdir; } ' EXIT

  asm=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE);
  samplename=$(basename $asm .fasta)
  b=$(basename $asm)
  echo "Sample name will be $samplename"

  # Speed this up by working on scratch
  cp -v $asm  $tmpdir/asm.fasta
  cp -v $db   $tmpdir/db.csv
  cp -v $refs $tmpdir/refs.fasta

  set -x
  EToKi.py MLSType -i $tmpdir/asm.fasta -r $tmpdir/refs.fasta -k $samplename -d $tmpdir/db.csv -o $tmpdir/out.fasta
  mv -v $tmpdir/out.fasta $outdir/$b

END_OF_SCRIPT

