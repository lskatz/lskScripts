#!/bin/bash
set -e
dir=$(realpath $(dirname $0));

lambda0=$($dir/../../scripts/Kendall.pl --alreadyrooted --lambda 0 $dir/input/kendall-colijn1.dnd $dir/input/kendall-colijn2.dnd | tail -n 1 | cut -f 4)
if [ "$lambda0" != "2.00" ]; then
  echo "For lambda0, got $lambda0 and not 2.00";
  exit 1;
fi;
lambda1=$($dir/../../scripts/Kendall.pl --alreadyrooted --lambda 1 $dir/input/kendall-colijn1.dnd $dir/input/kendall-colijn2.dnd | tail -n 1 | cut -f 4)
if [ "$lambda1" != "1.96" ]; then
  echo "For lambda1, got $lambda1 and not 1.96";
  exit 1;
fi;

echo "Success: $0";
exit 0;
