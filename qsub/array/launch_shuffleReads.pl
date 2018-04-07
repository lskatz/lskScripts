#!/bin/bash

set -e

if [ "$2" == "" ]; then
  echo "Only give _R1_ files to this script"
  echo "Usage:";
  exit 1
fi

outdir=$1
shift

TMP=$(mktemp --tmpdir='.' --directory SHUFFLE.XXXXXX)
CTRL_FILE="$TMP/fasta.txt"
echo "$@" | tr ' ' '\n' > $CTRL_FILE
mkdir -p $TMP/log

export PATH=$PATH:$HOME/bin/cg_pipeline/scripts

which run_assembly_shuffleReads.pl

qsub -q all.q -N shuffleReads -o $TMP/log -j y -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "outdir=$outdir" <<- "END_OF_SCRIPT"
  #!/bin/bash

  set -e
  
  R1=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  b=$(basename $R1);
  d=$(dirname  $R1);

  R2=${R1/_R1_/_R2_}

  echo "Shuffling on $(hostname) $R1 and $R2"

  ls $R1 $R2

  if [ -e "$outdir/$b.fastq.gz" ]; then
    echo "Found $outdir/$b.fastq.gz. Exiting.";
    exit 0;
  fi

  mkdir -p /dev/shm/$USER
  tmpdir=$(mktemp --tmpdir=/dev/shm/$USER --directory shuffleReads.XXXXXX);
  trap "{ rm -rf $tmpdir; }" EXIT

  shuffled="$tmpdir/$b.fastq.gz";
  
  run_assembly_shuffleReads.pl $R1 $R2 | gzip -9c > $shuffled

  mv -v $shuffled $outdir/

END_OF_SCRIPT

