#!/bin/bash 

# Runs mash sketch on a set of fasta files


# Read ARGV
FASTA=$@

set -e
set -u

if [ "$FASTA" == "" ]; then
  echo "Mash sketch files in place"
  echo "Usage: $0 *.fasta"
  exit 1;
fi

TMP=$(mktemp --tmpdir='.' --directory qsubMash.XXXXXXXX)
mkdir -p $TMP/log
echo "tmp dir is $TMP "

CTRL_FILE="$TMP/array.txt"
echo "$FASTA" | tr ' ' '\n'  > $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"

source /etc/profile.d/modules.sh
module purge
module load Mash
module load gcc/5.5
mash --version # ensure it loaded

qsub -q edlb.q -N mashSketch -o $TMP/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash -l

  set -e
  set -u
  set -x

  hostname
  # Reload modules to ensure things like LD_LIBRARY_PATH are re-added
  source /etc/profile.d/modules.sh || true
  module purge
  module load Mash
  module load gcc/5.5

  mash --version # ensure it loaded

  # Set up filenames
  fasta=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE);
  mkdir /scratch/$USER || true
  tmpdir=$(mktemp --tmpdir=/scratch/$USER --directory qsubMash.XXXXXXXX)
  trap ' { rm -rf $tmpdir; } ' EXIT

  tmpsketch="$tmpdir/out.msh"
  outsketch="$(dirname $fasta)/$(basename $fasta).msh"

  echo "mash will be run on $fasta => $tmpsketch => $outsketch"

  if [ -e "$outsketch" ]; then
    echo "$outsketch already exists. Exiting.";
    exit 1;
  fi

  mash sketch -o $tmpsketch $fasta

  mv -v $tmpsketch $outsketch
END_OF_SCRIPT

