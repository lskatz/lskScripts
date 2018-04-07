#!/bin/bash
set -e

NUMCPUS=1
DIR=$1

if [ "$DIR" == "" ]; then
  echo "Finds md5sum of all files in a whole directory, recursively"
  echo "Usage: $0 dir"
  exit 1
fi

# Make a temporary directory that will be removed upon exit
export TEMPDIR=$(mktemp --directory $(basename $0).XXXXXX --tmpdir)
function cleanup(){
  rm -rf $TEMPDIR
}
trap cleanup EXIT

# Find all files in the directory recursively,
# and md5sum them into TEMPDIR.
# Using the temporary directory helps with race conditions.
find $DIR -type f | xargs -P $NUMCPUS -n 1 sh -c '
  md5sum $0 > $TEMPDIR/$(basename $0).$$.md5
'
# Combine all md5sums to one file.
# Adding $$ to the filename will make it unique in the temp
# directory; the full file path appears in the md5 file
# contents and so it should be a stable sort.
sort $TEMPDIR/*.md5 | md5sum | cut -f 1 -d ' '
  
