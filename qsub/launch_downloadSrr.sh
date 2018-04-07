#!/bin/bash -l
#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o fastq-dump.log 
#$ -j y
#$ -N fastq-dump
#$ -q all.q

name=$1
OUT=$2

if [ "$OUT" == "" ]; then
  echo "Usage: $0 NAME"
  echo "  Downloads a set of fastq files using fastq-dump";
  exit 1;
fi;

module load sratoolkit/2.4.5-2
if [ $? -gt 0 ]; then echo "unable to load sratoolkit/2.4.5-2"; exit 1; fi;

mkdir -p /scratch/$USER/fastq-dump/$name

downloadSra.pl -t /scratch/$USER/fastq-dump/$name $name | gzip -c > $OUT
if [ $? -gt 0 ]; then echo "ERROR with fastq-dump $name"; exit 1; fi;

rm -rf /scratch/$USER/fastq-dump/$name

