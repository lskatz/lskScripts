#!/bin/bash

source /etc/profile.d/modules.sh
if [ $? -gt 0 ]; then 
  echo "ERROR: cannot load the modules system";
  exit 1;
fi

module purge

if [ "$3" == "" ]; then
  echo "Usage: $0 out.tsv in1.fasta in2.fasta [in3.fasta...]"
  echo "  Runs ANI on all vs all and places that output in out.tsv"
  exit
fi

OUT=$1
shift

if [ -e "$OUT" ]; then
  echo "ERROR: $OUT already exists"
  exit 1
fi

LOGDIR=$(mktemp --directory $(basename $0 .sh).XXXXXX)
CTRL_FILE="$LOGDIR/array.txt"
mkdir $LOGDIR/log
mkdir $LOGDIR/out
echo "log directory is $LOGDIR/log"

echo "$@" | tr ' ' '\n' > $CTRL_FILE

qsub -q edlb.q -q all.q -N FastANIarray -o $LOGDIR/out -e $LOGDIR/log -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash
  set -e

  QUERY=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)

  echo "ANI for query $QUERY" >&2
  hostname >&2

  fastANI -q $QUERY --rl $CTRL_FILE -o /dev/stdout

END_OF_SCRIPT

qsub -q all.q -N combine_FastANI -o $LOGDIR -j y -pe smp 1 -V -cwd -hold_jid FastANIarray \
  -v "LOGDIR=$LOGDIR" -v "OUT=$OUT" <<- "END_OF_SCRIPT"
  #!/bin/bash
  set -e

  sort -k1,2r $LOGDIR/out/*.o* | uniq > $OUT
END_OF_SCRIPT
