#!/bin/bash -l

#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o fastq_to_fasta.log -j y
#$ -N fastq_to_fasta
#$ -q all.q

in=$1
out=$2

if [ "$out" == "" ]; then
  echo "Usage: $0 in.fastq[.gz] out.fasta"
  exit 1;
fi;

extension="${in##*.}"

if [ "$extension" = "gz" ]; then
  gunzip -c "$in" | fastq_to_fasta -Q33 > "$out"
else
  fastq_to_fasta -Q33 < "$in" > "$out"
fi
if [ $? -gt 0 ]; then echo "ERROR with fastq_to_fasta"; exit 1; fi;
