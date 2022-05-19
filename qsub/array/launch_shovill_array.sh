#!/bin/bash -l

# Runs shovill on a set of Illumina reads


# Read ARGV
OUTDIR=$1; shift;
READS=$@

set -e

if [ "$READS" == "" ]; then
  echo "Assemble all reads in a directory from Illumina runs"
  echo "All files must be in order for paired-end to work"
  echo "Usage: $0 outdir *.fastq.gz"
  exit 1;
fi

if [ ! -d "$OUTDIR" ]; then
  mkdir "$OUTDIR"
fi;

tmpdir=$(mktemp --tmpdir='.' --directory qsubShovillSpades.XXXXXXXX)
mkdir -pv $tmpdir/log
#trap " rm -rf $tmpdir " EXIT
echo "tmp dir is $tmpdir "

# CTRL file will have per line:
#   filename  coverageLevel
CTRL_FILE="$tmpdir/array.txt"
# Put the reads one at a time into the CTRL_FILE but use paste to keep paired ends together
echo "$READS" | tr ' ' '\n' | paste - -  > $CTRL_FILE
#echo "DEBUG"; head $CTRL_FILE > $CTRL_FILE.tmp && mv $CTRL_FILE.tmp $CTRL_FILE
#echo "CTRL_FILE is $CTRL_FILE"

head -n 5 $CTRL_FILE
echo "This is what the top of the CTRL file looks like. "
echo "  Waiting 1 second in case you want to ctrl-c..."

# Might as well start loading the environment before sleeping
module purge
export PATH=$HOME/bin/shovill-v1.1.0/bin:$PATH
module load SPAdes/3.13.0 Skesa/2.3.0 megahit/1.1.2 velvet/1.2.10 lighter/1.1.1 flash/1.2.11 samtools/1.9 bwa/0.7.17 seqtk/1.3 pilon/1.22 trimmomatic/0.35 perl/5.16.1-MT kmc/3.0 java/jdk1.8.0_301
sleep 1

# See if we have all the right components
shovill --check

qsub -q all.q -q edlb.q -N ShovillSpades -o $tmpdir/log -j y -pe smp 4-6 -l h_vmem=72G -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "tmpdir=$tmpdir" -v "CTRL_FILE=$CTRL_FILE" -v "OUTDIR=$OUTDIR" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e
  set -u

  # bring back in modules that load LD_LIBRARY_PATH
  # because that variable is stripped away for security purposes
  module load gcc/4.9.3 pilon/1.22

  # Set up filenames
  fastq=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE);
  R1=$(echo "$fastq" | cut -f 1);
  R2=$(echo "$fastq" | cut -f 2);
  name=$(basename $R1 .fastq.gz)

  sampledir=$(mktemp --tmpdir=$tmpdir --directory $name.shovillSpades.XXXXXX);
  trap "rm -rf $sampledir" EXIT
  fastaOut="$OUTDIR/$name.shovillSpades.fasta"

  ls -lh $R1 $R2
  echo "Shovill with SPAdes will be run on $R1 $R2, on $(hostname)"
  echo "Fasta will be written to $fastaOut";

  if [ -e "$fastaOut" ]; then
    echo "$fastaOut already exists. Exiting.";
    exit 1;
  fi

  # Make sure the skesa lib is there before the full shovill check
  skesa --version 2>&1

  shovill --check

  shovill --R1 $R1 --R2 $R2 --outdir $sampledir --assembler spades --cpus $NSLOTS --force --tmpdir /scratch --ram 64
  cp -v $sampledir/contigs.fa $fastaOut
  # trap command will remove $sampledir

END_OF_SCRIPT

