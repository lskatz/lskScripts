#!/bin/bash

dir=$BATS_TEST_DIRNAME
export PATH=$dir/../../scripts:$PATH

@test "kendall" {

  lambda0=$($dir/../../scripts/Kendall.pl --alreadyrooted --lambda 0 $dir/input/kendall-colijn1.dnd $dir/input/kendall-colijn2.dnd | tail -n 1 | cut -f 4)
  [ "$lambda0" == "2.00" ]

  lambda1=$($dir/../../scripts/Kendall.pl --alreadyrooted --lambda 1 $dir/input/kendall-colijn1.dnd $dir/input/kendall-colijn2.dnd | tail -n 1 | cut -f 4)
  [ "$lambda1" == "1.96" ]
}
