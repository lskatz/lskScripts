#!/bin/bash

tree1=$1
tree2=$2

script=$(basename $0);
if [ "$tree1" == "" ]; then
  echo "Usage: $script tree1.dnd tree2.dnd"
  exit 1;
fi;

logmsg () {
  echo "$script: $@" >&2
}

# Check executables
for exe in treedist randTrees.pl; do
  which $exe >/dev/null 2>&1;
  if [ $? -gt 0 ]; then
    logmsg "ERROR: I could not find $exe in your path!";
    exit 1;
  fi;
done;

mkdir -p /tmp/$USER
tmpdir=$(mktemp --directory --tmpdir=/tmp/$USER Robinson-Foulds.XXXXXX)
if [ $? -gt 0 ]; then logmsg "ERROR making temporary directory under /tmp/$USER"; exit 1; fi;
logmsg "Temporary dir is $tmpdir";

cp $tree1 $tmpdir/ &&\
cp $tree2 $tmpdir/

t1="$tmpdir/$(basename $tree1)"
t2="$tmpdir/$(basename $tree2)"

# Create a list of trees to compare against:
#   The first tree is the tree to compare against.
#   The next trees are randomly made trees from 
#   the the tree to compare against's parameters
comparisonTrees="$tmpdir/compareAgainst.dnd";
cp $t2 $comparisonTrees
randTrees.pl $t1 --numTrees 100 >> $comparisonTrees
ln -s $t1 "$tmpdir/intree"
ln -s $comparisonTrees "$tmpdir/intree2"

# Must do treedist in the temp directory because it pollutes the
# local directory
pushd $tmpdir > /dev/null 2>&1

# Run treedist with the Kuhner and Felsenstein distance metric.
# If the option "D" were given, then it would be run with Robinson-Foulds
#echo -ne 'D\n2\nL\nS\nY\n' | treedist > "$tmpdir/treedist.log" 2> "$tmpdir/treedist.out"
echo -ne '2\nL\nS\nY\n' | treedist > "$tmpdir/treedist.log" 2> "$tmpdir/treedist.out"
if [ $? -gt 0 ]; then logmsg "ERROR with treedist program: $(cat $tmpdir/treedist.log)"; exit 1; fi;

# Find average and stdev
cat outfile | perl -MStatistics::Descriptive -MMath::Gauss=cdf,pdf -MList::Util=sum -lane '
  BEGIN{
    my @F=split(/\s+/,<>);
    $obs=$F[2];
  }
  next if($F[0] != 1);
  push(@difference,$F[2]);
  END{
    $stat=Statistics::Descriptive::Full->new();
    $stat->add_data(@difference);
    $num=@difference;
    $avg=$stat->mean;
    $stdev=$stat->standard_deviation;
    $var=$stdev**2;

    $Z=($obs - $avg)/$stdev;
    $p=cdf($obs,$avg,$stdev);
    #$p=cdf($Z);

    # scientific and floating point formatting
    $_=sprintf("%0.2e",$_) for($p);
    $_=sprintf("%0.2f",$_) for($obs,$avg,$stdev,$Z);

    # Print results
    print join("\t",qw(obs num avg stdev Z p));
    print join("\t",$obs, $num, $avg,$stdev,$Z,$p);
  }
'

rm -rf $tmpdir

