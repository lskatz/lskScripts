#!/bin/bash

set -e
set -x

# Find all unit tests under this directory
# and simply run them with -e and -x
find $(dirname $0)/unittests -maxdepth 1 -type f | \
  while read exe; do
    $exe
  done;

exit 0
