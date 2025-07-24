#!/bin/bash -l

set -e

if [[ -z $2 ]]; then
  echo "Runs quast on assemblies"
  echo "Usage: $0 outdir/ *.fasta"
  exit 0
fi

set -u

# get the output directory and then remove from argv
outdir=$1
shift

# Make temp folder that can hold log files and temporary files
TMP=$(mktemp --tmpdir='.' --directory quast.XXXXXXXX)
mkdir $TMP/log
echo "tmp dir is $TMP "
# Make a temp file that can hold an array of input files
CTRL_FILE="$TMP/array.txt"
echo $@ | tr ' ' '\n' | grep . > $CTRL_FILE

# Start off the job array
# -N is the job name
# -o $TMP -j y puts all log files into the temporary directory
# -V -cwd is to use the current environment in the current working directory
# -t indicates an array, but it needs a min to max in the next parameter
#     1-$(cat $CTRL_FILE | wc -l) translates to "1 to the number of gzip files"
# -v "CTRL_FILE=$CTRL_FILE" creates a variable to use inside of the here document.
qsub -pe smp 1 -N quast -o $TMP/log -j y -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "outdir=$outdir" <<- "END_OF_SCRIPT"
    # This is a "here document."  It gets submitted as though it were a 
    # separate file. The here document ends right before END_OF_SCRIPT

    set -e
    module load quast
    which quast.py

    file=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)

    TMP=$(mktemp --directory --suffix _qsub)
    echo "TMP DIR is $HOSTNAME:$TMP"
    trap " { rm -rf $TMP; } " EXIT
    tmpIn=$TMP/$(basename $file)
    uncompressed=$TMP/$(basename $file .gz)
    cp -v $file $tmpIn

    quast.py --threads 1 --output-dir $TMP/quast_out \
      --no-plots --no-html --no-html --no-icarus 

    mv -v $TMP/ $outdir/$(basename $file .fasta)_quast_out
END_OF_SCRIPT

