#!/bin/bash

export ref_tree=$1
export query_tree=$2

script=$(basename $0);
if [ "$ref_tree" == "" ]; then
  echo "Usage: $script ref.dnd query.dnd"
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
tmpdir=$(mktemp --directory --tmpdir=/tmp/$USER Kuhner-Felsenstein.XXXXXX)
if [ $? -gt 0 ]; then logmsg "ERROR making temporary directory under /tmp/$USER"; exit 1; fi;
logmsg "Temporary dir is $tmpdir";

cp $ref_tree $tmpdir/ &&\
cp $query_tree $tmpdir/

tmpRefTree="$tmpdir/$(basename $ref_tree)"
tmpQueryTree="$tmpdir/$(basename $query_tree)"

# Create a list of trees to compare against in $comparisonTrees:
#   The first tree is the reference tree
#   The next trees are randomly made trees from 
#   the query tree.
#   Therefore, the question being answered is, is the query
#   closer to the ref than random trees?
comparisonTrees="$tmpdir/compareAgainst.dnd";
cp $tmpRefTree $comparisonTrees
randTrees.pl $tmpQueryTree --numTrees 1000 >> $comparisonTrees
# We need the files to be named intree and intree2 because of inflexible treedist
ln -s $tmpQueryTree "$tmpdir/intree"
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

    # scientific and floating point formatting
    $_=sprintf("%0.2e",$_) for($p);
    $_=sprintf("%0.2f",$_) for($obs,$avg,$stdev,$Z);

    # Print results
    print join("\t",qw(ref_tree query_tree obs num avg stdev Z p));
    print join("\t",$ENV{ref_tree},$ENV{query_tree},$obs, $num, $avg,$stdev,$Z,$p);
  }
'

rm -rf $tmpdir

