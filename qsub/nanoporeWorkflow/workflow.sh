#!/bin/bash
#$ -o wtdbg2.log
#$ -j y
#$ -N wtdbg2
#$ -pe smp 1
#$ -V -cwd
#$ -l gpu=1
set -e

source /etc/profile.d/modules.sh
module purge

NSLOTS=${NSLOTS:=24}

OUTDIR=$1
FAST5DIR=$2
GENOMELENGTH=5000000 # TODO make this a parameter
LONGREADCOVERAGE=50  # How much coverage to target with long reads

set -u

thisDir=$(dirname $0);
thisScript=$(basename $0);

if [ "$FAST5DIR" == "" ]; then
    echo "Usage: $thisScript outdir fast5dir/"
    exit 1;
fi;

# Setup any debugging information
date
hostname

# Setup tempdir
tmpdir=$(mktemp -p . -d ONT-ASM.XXXXXX)
trap ' { echo "END - $(date)"; rm -rf $tmpdir; } ' EXIT
mkdir $tmpdir/log
echo "$0: temp dir is $tmpdir";

uuid1=$(uuidgen)
jobName1="basecall-$uuid1"
qsub -pe smp 1-$NSLOTS -N $jobName1 -cwd -o log/$jobName1.log -j y \
  $thisDir/01_basecall.sh $OUTDIR $FAST5DIR

# Now that it is demultiplexed, deal with each sample at a time.
for barcodeDir in $OUTDIR/demux/barcode[0-9]*; do

  # Prep the sample
  uuid2=$(uuidgen)
  jobName2="prepSample-$uuid2"
  qsub -hold_jid $jobName1 -pe smp 1-$NSLOTS -N $jobName2 -cwd -o log/$jobName2.log -j y \
    $thisDir/03_prepSample.sh $barcodeDir
  
  # Assemble the sample
  uuid3=$(uuidgen)
  jobName3="assemble-$uuid3"
  qsub -hold_jid $jobName2 -pe smp 1-$NSLOTS -N $jobName3 -cwd -o log/$jobName3.log -j y \
    $thisDir/05_assemble.sh $barcodeDir

  # Polish the sample
  uuid4=$(uuidgen)
  jobName4="assemble-$uuid4"
  qsub -hold_jid $jobName3 -pe smp 1-$NSLOTS -N $jobName4 -cwd -o log/$jobName4.log -j y \
    $thisDir/07_nanopolish.sh $barcodeDir $FAST5DIR

done
