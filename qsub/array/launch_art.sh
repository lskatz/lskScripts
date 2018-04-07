#!/bin/bash

# Launches ART whole genome sequence simulator on a set
# of fasta files
set -e

if [ "$4" == "" ]; then
  echo "Simulates reads from a set of fasta files."
  echo "Usage: $(basename $0) profileR1.txt profileR2.txt outdir/ file1.fasta file2.fasta..."
  exit 1
fi

profile1=$1;
profile2=$2
outdir=$3
shift;shift;shift;

# Check params
if [ ! -e $profile1 ]; then
  echo "ERROR: could not find $profile1"
  exit 1;
fi
if [ ! -e $profile2 ]; then
  echo "ERROR: could not find $profile2"
  exit 1;
fi
if [ -e $outdir ]; then
  echo "ERROR: outdir already exists: $outdir"
  exit 1;
fi
mkdir $outdir

export PATH=$PATH:~/bin/ART-MountRainier

TMP=$(mktemp --tmpdir='.' --directory ART.XXXXXX)
CTRL_FILE="$TMP/fasta.txt"
echo "$@" | tr ' ' '\n' > $CTRL_FILE

mkdir -p $TMP/log
qsub -q all.q -N ART -o $TMP/log -j y -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "outdir=$outdir" -v "profile1=$profile1" -v "profile2=$profile2" <<- "END_OF_SCRIPT"
  #!/bin/bash

  set -e
  
  fasta=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  b=$(basename $fasta | sed 's/\.[^.]*$//')

  if [ -e "$outdir/$b.1.fq.gz" ]; then
    echo "Found $outdir/$b.1.fq.gz. Exiting.";
    exit 0;
  fi

  mkdir -p /dev/shm/$USER
  tmpdir=$(mktemp --tmpdir=/dev/shm/$USER --directory ART.XXXXXX);
  trap "{ rm -rf $tmpdir; }" EXIT

  prefix="$tmpdir/$b." # will generate fastq with correct basename followed by '.1.fq.gz' or '.2.fq.gz'

  art_illumina -1 $profile1 -2 $profile2 -na -p -i $fasta -l 250 -f 40 -m 480 -s 120 -o $prefix
  gzip -v9 $prefix*

  mv -v $prefix*.fq.gz $outdir

END_OF_SCRIPT

