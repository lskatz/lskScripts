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

cp $tree1 $tree2 $tmpdir/

t1="$tmpdir/$(basename $tree1)"
t2="$tmpdir/$(basename $tree2)"

# Create a list of trees to compare against:
#   The first tree is the tree to compare against.
#   The next trees are randomly made trees from 
#   the the tree to compare against's parameters
comparisonTrees="$tmpdir/compareAgainst.dnd";
cp $t2 $comparisonTrees
randTrees.pl $t2 --numTrees 100 >> $comparisonTrees
ln -s $t1 "$tmpdir/intree"
ln -s $comparisonTrees "$tmpdir/intree2"

# Must do treedist in the temp directory because it pollutes the
# local directory
pushd $tmpdir > /dev/null 2>&1
echo -ne 'D\n2\nL\nS\nY\nY\n' | treedist > "$tmpdir/treedist.log" 2> "$tmpdir/treedist.out"
if [ $? -gt 0 ]; then logmsg "ERROR with treedist program: $(cat $tmpdir/treedist.log)"; exit 1; fi;

# Find average and stdev
total=$(cut -f 3 -d ' ' "$tmpdir/outfile")
cat outfile | perl -MList::Util=sum -lane '
  BEGIN{
    my @F=split(/\s+/,<>);
    $obs=$F[2];
  }
  next if($F[0] != 1);
  push(@difference,$F[2]);
  END{
    $num=@difference;
    $total=sum(@difference);
    $avg=($total/$num);
    
    $SS=0;
    for(@difference){
      $SS = $SS + ($_-$avg)**2;
    }
    $stdev=sqrt($SS/($num-1));
    $var=$stdev**2;

    $Z=($obs - $avg)/$stdev;
    $p=exp(-($obs - $avg)**2/(2 * $var))/sqrt(8 * atan2(1,1) * $var);

    $_=sprintf("%0.2f",$_) for($avg,$stdev,$Z,$p);
    print join("\t",qw(obs num avg stdev Z p));
    print join("\t",$obs, $num, $avg,$stdev,$Z,$p);
  }
'

rm -rf $tmpdir
