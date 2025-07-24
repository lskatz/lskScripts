#!/bin/bash -l
# taken from Seth's example: http://git.biotech.cdc.gov/snippets/3
# This is meant as an example to help me learn job arrays

# exit on error
set -e

# Print usage if no arguments
if [[ -z $1 ]]; then
  echo "Gzips fastq files with -9"
  echo "Usage: $0 *.fastq *.fastq.gz"
  exit 0
fi

# make an error with unset variables
set -u

# Make temp folder that can hold log files and temporary files
TMP=$(mktemp --tmpdir='.' --directory qsubtmp.XXXXXXXX)
mkdir $TMP/log
echo "tmp dir is $TMP "
# Make a temp file that can hold an array of input files
CTRL_FILE="$TMP/array.txt"
# put all command line arguments (fastq files) into the control file, one per line.
# Each SGE task below will have one specific line that it reads from.
echo $@ | tr ' ' '\n' | grep . > $CTRL_FILE

# Start off the job array
# -N is the job name
# -o $TMP -j y puts all log files into the temporary directory
# -V -cwd is to use the current environment in the current working directory
# -t indicates an array, but it needs a min to max in the next parameter
#     1-$(cat $CTRL_FILE | wc -l) translates to "1 to the number of gzip files"
# -v "CTRL_FILE=$CTRL_FILE" creates a variable to use inside of the here document.
qsub -N gzip-job -o $TMP/log -j y -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
    # This is a "here document."  It gets submitted as though it were a 
    # separate file. The here document ends right before END_OF_SCRIPT

    # the sed statement takes the line corresponding to the unique SGE task ID in the control file,
    # and the control file simply is a file of filenames, one per line.
    file=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)

    # the temp directory is used because it gives a unique space for this specific job
    # and many times it is made on a faster hard drive, making the analysis faster.
    TMP=$(mktemp --directory --suffix _qsub)
    echo "TMP DIR is $HOSTNAME:$TMP"
    # The trap statement cleans up the tmp drive at the end.
    # The -x is optional and lets you see the command(s) that run after it is called.
    trap " { set -x; rm -rf $TMP; } " EXIT

    # Start staging the gzip

    # tmpIn will be where the uncompressed file will be copied to temporarily
    tmpIn=$TMP/$(basename $file)
    uncompressed=$TMP/$(basename $file .gz)
    # copy the input file to our tmp space where it is in its own sandbox and will not collide with other SGE tasks
    cp -v $file $tmpIn

    if [[ "$tmpIn" =~ gz$ ]]; then
      gunzip -v $tmpIn && \
        gzip -9v $uncompressed && \
        mv -v $tmpIn $file
    else
      gzip -v9 $tmpIn && \
        mv -v $tmpIn.gz $file.gz && \
        rm -v $file
    fi
# end the HERE doc
END_OF_SCRIPT

