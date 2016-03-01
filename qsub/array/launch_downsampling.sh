#!/bin/bash -l

# downsamples a set of reads

if [ "$4" == "" ]; then
  echo "Usage: $0 oneX minCov maxCov 1.fastq.gz [2.fastq.gz ...]"
  echo "Reads will be deposited in cov\$i to cov\$j directories"
  echo "oneX must equal to the size of the reference assembly";
  echo "Can be fastq.gz or fastq files"
  exit 1;
fi

TMP=$(mktemp --tmpdir='.' --directory qsubDownsample.XXXXXXXX)
mkdir -p $TMP/log
echo "tmp dir is $TMP "

# Read ARGV
oneX=$1;
shift;
MIN=$1
MAX=$2
shift; shift;

# CTRL file will have per line:
#   filename  coverageLevel
CTRL_FILE="$TMP/array.txt"
(for cov in `seq $MIN $MAX`; do
  for j in "$@"; do
    echo "$j $cov";
  done;
done;) | grep . > $CTRL_FILE

qsub -q all.q -N downsample -o $TMP/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v oneX=$oneX -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash

  set -e

  # The global temporary directory
  tmpdir=/scratch

  echo "Downsampling script will be $(which run_assembly_removeDuplicateReads.pl)"

  # What coverage?  Directories?
  fastq=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk '{print $1}')
  b=$(basename $fastq);
  cov=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk '{print $2}')
  localDir=cov$cov/reads
  scratchDir=$tmpdir/$USER/cov$cov
  mkdir -p $scratchDir $localDir

  if [ ! -d $scratchDir ]; then
    echo "ERROR: could not make $scratchDir";
    exit 1;
  fi;

  # Final number of bp in the fastq file
  bp=$(($oneX * $cov));

  echo "Reading $fastq to a coverage of $cov ($bp bp)";

  tmpFastq=$scratchDir/$b;
  outFastq=$localDir/$b;
  run_assembly_removeDuplicateReads.pl $fastq --sizeto $bp --nobin | gzip -c > $tmpFastq 
  if [ $? -gt 0 ]; then
    echo "ERROR: $b" >> $localDir/ERROR
  fi
  mv -v $tmpFastq $outFastq

  # Clean up temporary directory if it's empty
  rmdir $scratchDir

END_OF_SCRIPT

