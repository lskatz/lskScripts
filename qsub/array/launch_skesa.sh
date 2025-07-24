#!/bin/bash -l

# Runs spades on a set of Illumina reads


# Read ARGV
OUTDIR=$1; shift;
READS=$@

set -e

if [ "$READS" == "" ]; then
  echo "Assemble all reads in a directory from Illumina runs"
  echo "All files must be in interleaved format";
  echo "Usage: $0 outdir *.fastq.gz"
  exit 1;
fi

if [ ! -d "$OUTDIR" ]; then
  mkdir "$OUTDIR"
fi;

TMP=$(mktemp --tmpdir='.' --directory qsubSkesa.XXXXXXXX)
mkdir -p $TMP/log
echo "tmp dir is $TMP "

# CTRL file will have per line:
#   filename  coverageLevel
CTRL_FILE="$TMP/array.txt"
echo "$READS" | tr ' ' '\n'  > $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"

module purge
module load Skesa
skesa --version

qsub -q all.q -N skesa -o $TMP/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "OUTDIR=$OUTDIR" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e

  module load Skesa

  module list
  skesa --version

  # Set up filenames
  fastq=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE);
  tmpdir=/scratch/$USER/skesa
  mkdir -p $tmpdir
  sampledir=$(mktemp --tmpdir=$tmpdir --directory skesa.XXXXXX);
  trap "rm -rf $sampledir" EXIT
  fastaOut="$sampledir/$(basename $fastq .fastq.gz).fasta"

  echo "skesa will be run on $fastq, on $(hostname)"
  echo "Temporary fasta will be written to $fastaOut";

  if [ -e "$fastaOut" ]; then
    echo "$fastaOut already exists. Exiting.";
    exit 1;
  fi

  /usr/bin/time -o $OUTDIR/time.$SGE_TASK_ID.tsv -f "$fastq\t%e" \
    skesa --cores $NSLOTS --fastq $fastq --gz --use_paired_ends > $fastaOut

  mv -v $fastaOut $OUTDIR
END_OF_SCRIPT

