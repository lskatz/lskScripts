#!/bin/bash 
#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o SRST.log -j y
#$ -N SRST2

source /etc/profile.d/modules.sh 
module load samtools/0.1.18;
module load bowtie2/2.1.0;

DB=$1
MLST_DEFS=$2
PREFIX=$3
OUTDIR=$4
INTERLEVED=$5

thisScript=`basename $0`;
if [ "$INTERLEVED" == "" ]; then
  echo "Usage: $thisScript MLSTDB.fasta MLST_DEFS.txt OUTPREFIX OUTDIR/ INTERLEVED.fastq[.gz]"
  echo "OUTPREFIX cannot include a prefix directory which is why you should specify OUTDIR"
  exit 1;
fi

NSLOTS=${NSLOTS-1}
SCRIPT=/scicomp/home/gzu2/bin/srst2/scripts/srst2.py

b=$(basename $INTERLEVED);
READ1="TMP/$b.read_1.fastq"
READ2="TMP/$b.read_2.fastq"
run_assembly_shuffleReads.pl -d $INTERLEVED > $READ1 2> $READ2
if [ $? -gt 0 ]; then exit 1; fi;

# Low compression to make it compatible with srst2.
# These files will be deleted later anyway
gzip -f -v -1 $READ1 $READ2
if [ $? -gt 0 ]; then exit 1; fi;

READ1="$READ1.gz"
READ2="$READ2.gz"


$SCRIPT --input_pe $READ1 $READ2 --mlst_delimiter "_" --output $PREFIX --log --mlst_db $DB --mlst_definitions $MLST_DEFS
if [ $? -gt 0 ]; then exit 1; fi;
mv -v ${PREFIX}*__results.txt "$OUTDIR/"
if [ $? -gt 0 ]; then exit 1; fi;

# remove temp files
rm -v $READ1 $READ2
rm -vf ${PREFIX}__*.bam* ${PREFIX}__*.sam* ${PREFIX}__*results.txt ${PREFIX}__*.pileup ${PREFIX}*.log
