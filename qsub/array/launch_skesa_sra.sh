#!/bin/bash -l

# Runs spades on a set of Illumina reads


# Read ARGV
OUTDIR=$1; shift;
READS=$@

set -e

if [ "$READS" == "" ]; then
  echo "Assemble all sra runs in a file"
  echo "Usage: $0 outdir allsra_acc.txt"
  echo "  where the sra accessions are separated by any whitespace"
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
cat $READS | perl -lane '
  chomp;
  for my $f (split(/\s+/, $_)) {
    print "$f";
  }
' > $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"

# Check to make sure skesa will work in the heredoc before getting there
module purge
module load Skesa
skesa --version
module purge

qsub -q all.q -N skesa -o $TMP/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "OUTDIR=$OUTDIR" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e

  module purge
  module load Skesa
  module list
  skesa --version

  # Set up filenames
  run_acc=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE);
  tmpdir=/scratch/$USER/skesa
  mkdir -pv $tmpdir
  sampledir=$(mktemp --tmpdir=$tmpdir --directory skesa.XXXXXX);
  trap "rm -rf $sampledir" EXIT
  fastaOut="$sampledir/$run_acc.skesa.fasta"
  finalOut="$OUTDIR/$run_acc.skesa.fasta"

  echo "skesa will be run on $run_acc, on $(hostname)"
  echo "Temporary fasta will be written to $fastaOut";

  if [ -e "$fastaOut" ]; then
    echo "$fastaOut already exists. Exiting.";
    exit 2;
  fi

  if [ -e "$finalOut" ]; then
    echo "$finalOut already exists. Exiting.";
    exit 3;
  fi

  skesa --cores $NSLOTS --sra_run $run_acc > $fastaOut

  mv -v $fastaOut $finalOut
END_OF_SCRIPT

