#!/bin/bash -l

# Runs any set of reads through Lyve-SET in a cluster-friendly way.
# Each reads directory will be a distinct project.
# Author: Lee Katz
# Usage: bash $0 reference.fasta readsdir1 readsdir2 [... readsdirX]

if [ "$2" == "" ]; then
  echo "Usage: $0 ref/ dir [dir2 ... ]"
  echo "  The reference directory can just contain reference.fasta or can have other Lyve-SET reference directory files."
  echo "  Each directory will be searched for shuffled reads files matching *.f*q.gz"
  exit 1;
fi

TMP=$(mktemp --tmpdir='.' --directory qsubLyveSET.XXXXXXXX)
echo "tmp dir is $TMP "

REF=$1; shift; # get the reference genome and remove it from ARGV
CTRL_FILE="$TMP/array.txt"
echo "$@" | tr ' ' '\n' | grep . > $CTRL_FILE

mkdir -p $TMP/log
qsub -q all.q -N LyveSetShuffled -o $TMP/log -j y -pe smp 3-4 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "REF=$REF" <<- "END_OF_SCRIPT"
  #!/bin/bash

  set -e

  base_dir=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  echo "Working on $base_dir"
  mkdir -p /scratch/$USER
  if [ -e "$base_dir/Lyve-SET" ]; then
    echo "Found $base_dir/Lyve-SET! Will not continue.";
    exit 0;
  fi;
  scratch_out=$(mktemp --tmpdir="/scratch/$USER" --directory Lyve-SET.XXXXXX)
  rm -rfv $scratch_out

  # Run Lyve-SET on the scratch drive without qsub inception.
  set_manage.pl --create $scratch_out
  rmdir $scratch_out/reference;
  cp -r $REF $scratch_out/reference
  ln -sv $(find $(realpath $base_dir) -name '*.f*q.gz') $scratch_out/reads/
  if [ $? -gt 0 ]; then exit 1; fi;
  launch_set.pl --noqsub --numcpus $NSLOTS --read_cleaner CGP --mask-phages --mask-cliffs $scratch_out

  rm -rvf $scratch_out/{reads,bam,tmp}/* # no need to take up all this space
  mv -v $scratch_out $base_dir/Lyve-SET
END_OF_SCRIPT

