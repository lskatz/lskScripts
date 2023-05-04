#!/bin/bash

dir=$BATS_TEST_DIRNAME
export PATH=$dir/../../scripts:$PATH

@test "randTrees" {
  tree=$($dir/../../scripts/randTrees.pl --numTrees 1 $dir/input/kendall-colijn1.dnd)
  bytes=$(wc -c <<< $tree)
  [ "$bytes" -ge 140 ]

  [ "$bytes" -le 150 ]

}
