#!/bin/bash -l

# Modified from a script by Taylor Griswold <ycj5@cdc.gov>

# Read ARGV
OUTDIR=$1
DOWNLOAD_LIST=$2
SLOTS_PER_JOB=1 # manually change this as needed

if [ "DOWNLOAD_LIST" == "" ]; then
  echo "Executes the sratoolkit module into an array batch job."
  echo "Text file has white-space delimited SRA run IDs"
  echo "Usage: $0 outdir run_ids.txt"
  exit 1;
fi

if [ ! -d $OUTDIR ]; then
  echo "ERROR: $OUTDIR is not a directory";
  exit 1;
fi
if [ ! -e "$DOWNLOAD_LIST" ]; then
  echo "ERROR: $DOWNLOAD_LIST could not be found";
  exit 1;
fi

TMP=$(mktemp --tmpdir='.' --directory qsubFastqDump.XXXXXXXX)
mkdir -p $TMP/log
echo "tmp dir is $TMP "

# CTRL file will have one SRA run ID per line
CTRL_FILE="$TMP/array.txt"
cat $DOWNLOAD_LIST | perl -lane '
  for my $sra(@F){
    print $sra;
  }
' > $CTRL_FILE
echo "CTRL_FILE is $CTRL_FILE"


qsub -N FastqDump -q edlb.q -o $TMP/log -j y -pe smp $SLOTS_PER_JOB -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "OUTDIR=$OUTDIR" -v "CTRL_FILE=$CTRL_FILE","DOWNLOAD_LIST=$DOWNLOAD_LIST" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e
  source /etc/profile.d/modules.sh
  module purge
  module load sratoolkit/2.9.1
  
  which fastq-dump

  # Set up filenames
  SRRID=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE | awk '{print $1}')
  READ1=$(echo "$SRRID"_1.fastq.gz)
  READ2=$(echo "$SRRID"_2.fastq.gz)
  tmpdir=/scratch/$USER
  mkdir -p $tmpdir
  outdir=$(mktemp --tmpdir=$tmpdir --directory FastqDump.XXXXXX);
  trap "rm -rf $outdir" EXIT
  newdir="$OUTDIR/$SRRID.fastq-dump"

  if [ -e "$newdir" ]; then
    echo "ERROR: found pre-existing dir $newdir. Will not continue.";
    exit 1
  fi

  echo "fastq-dump will be run on $SRRID under $(hostname)"
  echo "Working directory is $outdir"
  echo -e "Final directory will be $newdir.\n";

  fastq-dump --accession $SRRID --outdir $outdir --defline-seq '@$ac.$si/$ri' --defline-qual '+' --split-files --skip-technical --dumpbase --clip --gzip
  
  mv -v $outdir $newdir
  rm /scicomp/home/ycj5/ncbi/public/sra/"$SRRID".sra
  
END_OF_SCRIPT

