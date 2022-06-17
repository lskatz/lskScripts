#!/bin/bash

set -e;

FASTQLIST="$1"
ref="$2"

ref=$(realpath $ref)
echo "Testing the path to the reference genome"
ls $ref

TMP=$(mktemp --tmpdir='.' --directory $(basename $0 .sh).tmp.XXXXXX)
echo "tmp dir is $TMP"
mkdir $TMP/samples
mkdir $TMP/log

CTRL_FILE="$TMP/array.txt"
grep '_1.*fastq.gz' $FASTQLIST > $CTRL_FILE

qsub -N genotypi -o $TMP/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "TMP=$TMP" -v "ref=$ref" <<- "END_OF_SCRIPT"
#/bin/bash -l
set -e
set -x

  R1=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  R2=${R1/_1/_2};
  out="out/$(basename $R1 .fastq.gz).tsv";
  out=$(realpath $out)
  if [ -e "$out" ]; then
    #echo "FOUND $out";
    exit 0
  fi


  set +x
  module load samtools/1.10 bcftools/1.10.2 bowtie2/2.3.5.1
  set -x

  export PATH=$PATH:$HOME/bin/genotyphi


  # I think that this script pollutes the working directory
  # but I might be wrong.  Still though, just change
  # the working directory to something in TMP to avoid problems.
  workingDir="$TMP/samples/$(basename $R1)"
  mkdir -pv $workingDir
  cd $workingDir
  bowtie2 -p $NSLOTS -x $ref -1 $R1 $R2 | samtools view -bh - > unsorted.bam
  samtools sort -o sorted.bam unsorted.bam
  samtools index sorted.bam
  genotyphi.py --mode bam --bam sorted.bam --ref $ref --ref_id AL513382.1 --output $(basename $out)

  # just in case I needed to see files before something broke
  ls -lht
  
  mv *.tsv $out

END_OF_SCRIPT

exit

#### original script

grep '_1.*fastq.gz' ../fofn.txt |\
  while read R1; do
    R2=${R1/_1/_2};
    bam="bam/$(basename $R1 .fastq.gz).bam";
    out="out/$(basename $R1 .fastq.gz).tsv";
    if [ -e "$out" ]; then
      echo "FOUND $out";
      continue;
    fi;

    bowtie2 -p 24 -x CT18.fasta -1 $R1 $R2 2> "$out.log" |\
      samtools view -bh - > unsorted.bam;
    samtools sort -o $bam unsorted.bam;
    python genotyphi/genotyphi.py --mode bam --bam $bam --ref CT18.fasta --ref_id AL513382.1 --output $(basename $out) >> "$out.log" 2>&1;
    localout=$(\ls *$(basename "$out"));
    mv -v "$localout" "$out";
  done

