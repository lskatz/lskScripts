#!/bin/bash -l
#$ -S /bin/bash
#$ -pe smp 4-16
#$ -cwd -V
#$ -o spades.log
#$ -j y
#$ -N SPAdes3.1.0
#$ -q all.q

reads=$1
out=$2
fastaOut=$3

if [ "$out" == "" ]; then
  echo "Usage: $0 reads.fastq.gz output/ [out.fasta]"
  echo "  if out.fasta is given, then the output directory will be removed and the scaffolds.fasta file will be saved"
  exit 1;
fi;

module load SPAdes/3.1.0
if [ $? -gt 0 ]; then echo "unable to load spades 3.1.0"; exit 1; fi;

NSLOTS=${NSLOTS:=1}
spades.py --12 $reads --careful -o $out -t $NSLOTS 
if [ $? -gt 0 ]; then echo "problem with spades 3.1.0"; exit 1; fi;

# Assembly metrics. Don't die if this script dies.  It's not worth it.
echo "# CG-Pipeline metrics" > $out/run_assembly_metrics.txt
run_assembly_metrics.pl $out/scaffolds.fasta >> $out/run_assembly_metrics.txt

if [ "$fastaOut" != "" ]; then
  cp -v "$out/scaffolds.fasta" $fastaOut
  if [ $? -gt 0 ]; then echo "problem with copying $out/scaffolds.fasta => $fastaOut"; exit 1; fi;
  rm -rf "$out";
  if [ $? -gt 0 ]; then echo "problem with removing the directory $out"; exit 1; fi;
fi

