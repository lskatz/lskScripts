#!/bin/bash -l

# Runs any TreeToReads projects in a cluster-friendly method
# Author: Lee Katz
# Usage: bash launch_TTR.sh project1 project2 [... projectX]
#   where each project has its own TTR.cfg file and associated TTR files.

TMP=$(mktemp --tmpdir='.' --directory qsubTTR.XXXXXXXX)
echo "tmp dir is $TMP "

CTRL_FILE="$TMP/array.txt"
echo "$@" | tr ' ' '\n' > $CTRL_FILE

mkdir -p $TMP/log
qsub -q all.q -N TTR -o $TMP/log -j y -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  export PATH=$PATH:~/bin/TreeToReads:~/bin/ART
  module unload perl/5.16.1-MT
  export PERL5LIB=""

  base_dir=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  echo "Working on $base_dir"
  scratch_out="/scratch/gzu2/TTR/$base_dir"
  rm -rfv $scratch_out $base_dir/out
  mkdir -p $(dirname $scratch_out)
  cd $base_dir
  sed -i.bak "s|output_dir.*|output_dir = $scratch_out|" TTR.cfg
  treetoreads.py TTR.cfg
  cd -
  mv -v $scratch_out $base_dir/out
END_OF_SCRIPT

