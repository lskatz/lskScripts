#!/bin/bash

dir=$BATS_TEST_DIRNAME
export PATH=$dir/../../scripts:$PATH

@test "avgstdev" {
  total="Total: 66"
  average="Average: 6.00 +/- 3.32"
  median="Median: 6.00 [3.50,8.50] [1.00-11.00]"
  mad="MAD: 3.00"

  observed=$(seq 1 11 | $dir/../../scripts/avgstdev.pl)

  [ "$(grep Total <<< "$observed")" == "$total" ]

  [ "$(grep Average <<< "$observed")" == "$average" ]

  [ "$(grep Median <<< "$observed")" == "$median" ]

  [ "$(grep MAD <<< "$observed")" == "$mad" ]

}
