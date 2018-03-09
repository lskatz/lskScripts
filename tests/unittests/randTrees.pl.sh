#!/bin/bash

set -e

tree=$(./randTrees.pl --numTrees 1 tests/unittests/input/kendall-colijn1.dnd)
bytes=$(wc -c <<< $tree)
if [ "$bytes" -lt 140 ]; then
  echo "Size of the tree is less than 140 bytes, weird. The tree was:"
  echo "  $tree"
  exit 1
fi
if [ "$bytes" -gt 150 ]; then
  echo "Size of the tree is greater than 150 bytes, weird. The tree was:"
  echo "  $tree"
  exit 1
fi

echo "randTrees.pl test passed!"
exit 0
