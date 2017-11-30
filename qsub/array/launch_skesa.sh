#!/bin/bash -l

# Runs spades on a set of Illumina reads

set -e

# Read ARGV
OUTDIR=$1; shift;
READS=$@

if [ "$READS" == "" ]; then
  echo "Assemble all reads in a directory from Illumina runs"
  echo "All files must be in interleaved format";
  echo "Usage: $0 outdir *.fastq.gz"
  exit 1;
fi

if [ ! -d "$OUTDIR" ]; then
  echo "ERROR: could not find out dir $OUTDIR";
  exit 1;
fi;

TMP=$(mktemp --tmpdir='.' --directory qsubSkesa.XXXXXXXX)
mkdir -p $TMP/log
echo "tmp dir is $TMP "

# CTRL file will have per line:
#   filename  coverageLevel
CTRL_FILE="$TMP/array.txt"
echo "$READS" | tr ' ' '\n' | grep "_R1_" > $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"

module purge
module load Skesa/2.0_2
#module unload gcc
#module load gcc/4.9.3
skesa --version

qsub -q all.q -N skesa -o $TMP/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "OUTDIR=$OUTDIR" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e

  # LD_LIBRARY_PATH is stripped by qsub and needs to be readded
  # The following is how it appeared in my own environment.
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/x86_64/gcc/5.4/lib64:/apps/x86_64/gmp/6.1.0/lib:/apps/x86_64/mpfr/3.1.3/lib:/apps/x86_64/mpc/1.0.3/lib:/apps/x86_64/isl/0.18/lib
  echo "WARNING: setting LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

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

  skesa --cores $NSLOTS --fastq $fastq --gz --use_paired_ends > $fastaOut

  mv -v $fastaOut $OUTDIR
END_OF_SCRIPT

