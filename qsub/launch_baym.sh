#!/bin/bash
# Runs a shovill assembly.  Run with no options for usage.
# Author: Lee Katz <lkatz@cdc.gov>
# Workflow taken from Dorian Feistel

#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o baym.log
#$ -j y
#$ -N baym
#$ -q all.q

outdir=$1; shift
reads=$@
NSLOTS=${NSLOTS:=1}

source /etc/profile.d/modules.sh
scriptname=$(basename $0);

if [ "$reads" == "" ]; then
  echo "Usage: $scriptname outdir *_1.fastq[.gz]"
  echo "  R2 reads will be detected automatically based on matchiing filenames"
  exit 0;
fi;

set -e
set -u

module purge
module load seqtk/1.3 trimmomatic/0.39 bwa/0.7.17 java/latest

# die if some crucial program is not present
# And side bonus of knowing where each executable is loading from
which seqtk
which bwa
which java
which kallisto
which 1output_abundances.py
which trim_reads_nwss.sh

reference_db=/scicomp/groups-pure/Projects/NWSS_SequencingData/apps/wastewater_analysis/NWSS_PIPELINE/03.reference_set_29SEP2021-2/sequences.kallisto_idx
scripts=/scicomp/groups-pure/Projects/NWSS_SequencingData/apps/wastewater_analysis/NWSS_PIPELINE/scripts

tmpdir=$(mktemp --directory --suffix=$(basename $0));
trap ' { rm -rf $tmpdir; } ' EXIT

mkdir -v $outdir || echo "WARNING: outdir already exists: $outdir"

for R1 in $reads; do
  name=$(basename $R1 .fastq.gz | perl -plane 's/_\d|\.fastq|\.gz//g');
  sampledir="$tmpdir/$name"
  mkdir -v $sampledir

  # Get the file extension, if it's .gz
  ext=${R1##*.}

  # This file has to be local
  ln -sv $(which 1output_abundances.py) $sampledir/

  # Get fastq files into our local tmp folder
  R2=${R1/_1.f/_2.f}
  # decompress or simply cat R1/R2 with zcat
  local_R1="$sampledir/$(basename $R1 .gz)"
  local_R2="$sampledir/$(basename $R2 .gz)"

  if [[ "$ext" == "gz" ]]; then
    zcat $R1 > $local_R1
    zcat $R2 > $local_R2
  else
    cp -vL $R1 $R2 $sampledir/
  fi

  # WARNING
  # Some commands require that you are in the sample directory
  # And so I will indent to indicate CWD
  cd $sampledir

    # Should just yield one sample name
    # into sample_ID
    local_sample_ID="sample_IDs.txt"
    \ls -f1 *.fastq | sed 's/_[12].fastq.*//' | sort | uniq > $local_sample_ID

    trim_dir="trimmed.out"
    kallisto="kallisto.out"

    echo "SAMPLE(S): $(cat $local_sample_ID | tr '\n' ',')"

    bash trim_reads_nwss.sh . $trim_dir $local_sample_ID

    echo

    bash run_kallisto_WG.sh $local_sample_ID $trim_dir/ivar $kallisto

  # Moving back out of the sample directory
  cd -

  cp -rv $sampledir/*kallisto.out $outdir/

  #run_kallisto_SPIKE.sh

  # Directly remove the directory since we're done with it now
  rm -rvf $sampledir
done

