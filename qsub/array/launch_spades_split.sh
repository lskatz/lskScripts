#!/bin/bash -l

# Runs spades on a set of Illumina reads

# Read ARGV
READS=$@

if [ "$READS" == "" ]; then
  echo "Assemble all reads in a directory from Illumina runs"
  echo "All files must be in the format *_R1_*.fastq.gz to ensure pairs stay together"
  echo "Usage: $0 *.fastq.gz"
  exit 1;
fi

TMP=$(mktemp --tmpdir='.' --directory qsubSPAdes.XXXXXXXX)
mkdir -p $TMP/log
echo "tmp dir is $TMP "

# CTRL file will have per line:
#   filename  coverageLevel
CTRL_FILE="$TMP/array.txt"
echo "$READS" | tr ' ' '\n' | grep "_R1_" > $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"

qsub -q all.q -N spades -o $TMP/log -j y -pe smp 1-99 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e
  module load SPAdes

  which spades.py

  # Set up filenames
  R1=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk '{print $1}')
  R2=${R1/_R1_/_R2_}
  tmpdir=/scratch/$USER
  mkdir -p $tmpdir
  outdir=$(mktemp --tmpdir=$tmpdir --directory SPAdes.XXXXXX);
  trap "rm -rf $outdir" EXIT
  newdir=$(dirname $R1)/$(basename $R1 .fastq.gz).spades

  if [ -e "$newdir" ]; then
    echo "ERROR: found pre-existing dir $newdir. Will not continue.";
    exit 1
  fi

  echo "spades will be run on $(hostname)"
  echo "R1/R2:"
  echo "  $R1"
  echo "  $R2"
  echo "Working directory is $outdir"
  echo "Final directory will be $newdir";

  spades.py -t $NSLOTS -1 $R1 -2 $R2 -o $outdir --careful

  mv -v $outdir $newdir
END_OF_SCRIPT

