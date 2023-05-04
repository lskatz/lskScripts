#!/bin/bash

set -e

dir=$(realpath $(dirname $0));
export PATH=$dir/../../scripts:$PATH

total="Total: 66"
average="Average: 6.00 +/- 3.32"
median="Median: 6.00 [3.50,8.50] [1.00-11.00]"
mad="MAD: 3.00"

observed=$(seq 1 11 | $dir/../../scripts/avgstdev.pl)

if [ "$(grep Total <<< "$observed")" != "$total" ]; then
  echo "Total needed to be $total"
  echo "Observed $observed"
  exit 1
fi

if [ "$(grep Average <<< "$observed")" != "$average" ]; then
  echo "Average needed to be $average"
  echo "Observed $observed"
  exit 1
fi
if [ "$(grep Median <<< "$observed")" != "$median" ]; then
  echo "Median needed to be $median"
  echo "Observed $observed"
  exit 1
fi
if [ "$(grep MAD <<< "$observed")" != "$mad" ]; then
  echo "MAD needed to be $mad"
  echo "Observed $observed"
  exit 1
fi

echo "avgstdev test passed!"
exit 0
