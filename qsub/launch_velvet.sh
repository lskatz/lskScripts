#!/bin/bash
#$ -pe smp 16
#$ -cwd -V
#$ -o velvet.log -j y
#$ -N velvet

module ()
{
  eval `/usr/bin/modulecmd bash $*`
}

module load velvet/1.2.10;
if [ $? -gt 0 ]; then echo "unable to load velvet/1.2.10"; exit 1; fi;

reads=$1
out=$2

# number of cpus is either set by SGE, or is just 1
#NSLOTS=${NSLOTS:=1}
NSLOTS=${NSLOTS:=1}
echo $NSLOTS

if [ "$out" == "" ]; then
  echo "Usage: $0 reads.fastq.gz output/"
  exit 1;
fi;

command="$(which perl) $(which VelvetOptimiser.pl) -s 55 -e 99 -d $out -p $out -t $NSLOTS -f '-fastq.gz -shortPaired $reads'"
eval $command
if [ $? -gt 0 ]; then echo "problem with VelvetOptimiser"; echo $command; exit 1; fi;
