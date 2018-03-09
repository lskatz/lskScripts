#!/bin/bash

# Transform a TTR simulation to a MiSeq run

if [ "$2" == "" ]; then
  echo "Transform a TTR directory to a MiSeq run"
  echo "Usage: $0 TTR out.miseq"
  exit 1;
fi

IN=$1
OUT=$2

mkdir -p /tmp/$USER
tmpdir=$(mktemp --directory --tmpdir=/tmp/$USER ttrToMiSeq.XXXXXX)


# Sample sheet
RUNNAME=$(basename IN)
READLENGTH=$(grep read_length $IN/TTR.cfg | grep -o [0-9]*)
if [ $? -gt 0 ]; then echo "ERROR reading $IN/TTR.cfg"; exit 1; fi;
DATE=$(date +'%m/%d/%Y')

CSV="$tmpdir/SampleSheet.csv"
echo -ne "[Header] 
IEMFileVersion,4
Investigator Name,TreeToReads
Experiment Name,$RUNNAME
Date,$DATE
Workflow,GenerateFASTQ
Application,FASTQ Only
Assay,Nextera XT
Description,$RUNNAME
Chemistry,Amplicon
  
[Reads] 
$READLENGTH
$READLENGTH

[Settings]  
ReverseComplement,0
Adapter ATCGATCGATCG

[Data]
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
" > $CSV

SAMPLE=$(ls $IN/fastq | xargs -n 1 basename)
for i in $SAMPLE; do
  echo -e "$i,,$RUNNAME,A01,N801,ATCGAAA,S801,ATCGAAA,1," >> $CSV
done

# Fastq files
FASTQDIR="$tmpdir/Data/Intensities/BaseCalls"
mkdir -p $FASTQDIR

for i in $SAMPLE; do
  set -e
  cp -v $IN/fastq/$i/*_1.fq.gz $FASTQDIR/${i}_S1_L001_R1_001.fastq.gz
  cp -v $IN/fastq/$i/*_2.fq.gz $FASTQDIR/${i}_S1_L001_R2_001.fastq.gz
  set +e
done;

echo "$tmpdir -> $OUT"
mv $tmpdir $OUT

