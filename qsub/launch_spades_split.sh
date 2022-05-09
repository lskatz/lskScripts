#!/bin/bash -l
# Runs a spades assembly.  Run with no options for usage.
# Author: Lee Katz <lkatz@cdc.gov>
#Example:  (for f in *R1_001.fastq.gz; do b=`basename $f _R1.fastq.gz`; r=`sed 's/R1/R2/' <<< $f`; qsub -N spades$b -o ./assemblies/log/b.spades.log ~/bin/launch_SPAdes_v3.11.0.sh $f $r ./assemblies/$b.spades3.11 ./assemblies/$b; done;)

#$ -S /bin/bash
#$ -pe smp 4-16
#$ -cwd -V
#$ -o spades.log
#$ -j y
#$ -N SPAdes3.15.3

forward=$1
reverse=$2
out=$3
fastaOut=$4

source /etc/profile.d/modules.sh
scriptname=$(basename $0);

if [ "$out" == "" ]; then
  echo "Usage: $scriptname reads.1.fastq.gz reads.2.fastq.gz output/ [out.fasta]"
  echo "  if out.fasta is given, then the output directory will be removed and the scaffolds.fasta file will be saved"
  exit 1;
fi;

module load SPAdes/3.15.3
if [ $? -gt 0 ]; then echo "unable to load spades 3.15.3"; exit 1; fi;

NSLOTS=${NSLOTS:=1}

COMMAND="spades.py -1 $forward -2 $reverse --careful -o $out -t $NSLOTS"
echo "$scriptname: $COMMAND"
$COMMAND
if [ $? -gt 0 ]; then echo "problem with spades 3.11.0"; exit 1; fi;

if [ "$fastaOut" != "" ]; then
  cp -v "$out/scaffolds.fasta" $fastaOut
  if [ $? -gt 0 ]; then echo "problem with copying $out/scaffolds.fasta => $fastaOut"; exit 1; fi;
  rm -rf "$out";
  if [ $? -gt 0 ]; then echo "problem with removing the directory $out"; exit 1; fi;
fi
