#!/bin/bash -l

#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o mergeFastaReads.log
#$ -j y
#$ -N mergeFastaReads
#$ -q all.q

reads=$1
out=$2

if [ "$out" == "" ]; then
  echo "Usage: $0 reads.fasta merged.fasta"
  exit 1;
fi;

merge_fasta_reads "$reads" > "$out.tmp" && mv "$out.tmp" "$out"
if [ $? -gt 0 ]; then echo "ERROR with merge_fasta_reads"; exit 1; fi;
