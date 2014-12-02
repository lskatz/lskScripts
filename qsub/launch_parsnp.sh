#!/bin/bash -l

#$ -S /bin/bash
#$ -pe smp 4-16
#$ -cwd -V
#$ -o parsnp.log
#$ -j y
#$ -N parsnp
#$ -q all.q

module load harvest/1.0.1
if [ $? -gt 0 ]; then echo "WARNING: couldn't load harvest module"; fi;

refGbk=$1
asmDir=$2
out=$3
script=$(basename $0)

if [ "$out" == "" ]; then
  echo "Usage: $script reference.gbk asmDir outDir"
  exit 1;
fi;

NSLOTS=${NSLOTS:=12}

c="parsnp -a 13 -c -R 1 -g $refGbk -d $asmDir -p $NSLOTS -o $out"
$c # run the command
if [ $? -gt 0 ]; then echo -e "ERROR with parsnp\n  $c"; exit 1; fi;
