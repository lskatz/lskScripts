#!/bin/bash -l

# Runs any set of reads through Snp-Pipeline in a cluster-friendly way.
# Each reads directory will be a distinct project.
# Author: Lee Katz
# Usage: bash $0 reference.fasta readsdir1 readsdir2 [... readsdirX]

TMP=$(mktemp --tmpdir='.' --directory qsubSnp-Pipeline.XXXXXXXX)
echo "tmp dir is $TMP "

REF=$1; shift; # get the reference genome and remove it from ARGV
CTRL_FILE="$TMP/array.txt"
echo "$@" | tr ' ' '\n' | grep . > $CTRL_FILE

mkdir -p $TMP/log
# Has to have an exclusive node because it is a greedy script
qsub -q all.q -N snp-pipeline -o $TMP/log -j y -pe smp 12,16 -hard -l exclusive=true -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" -v "REF=$REF" <<- "END_OF_SCRIPT"
  #!/bin/bash

  set -e

  module load bowtie2/2.1.0
  module load varscan/2.3.7
  export PATH=~/.local/bin:$PATH # make sure snp-pipeline is prioritized
  export CLASSPATH=/apps/x86_64/varscan/bin/VarScan.v2.3.7.jar # varscan

  base_dir=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  echo "Working on $base_dir"
  mkdir -p /scratch/$USER
  scratch_out=$(mktemp --tmpdir="/scratch/$USER" --directory snp-pipeline.XXXXXX)
  export scratch_out
  mkdir -p $scratch_out/samples

  # Make the config file and put it into scratch_out (snppipeline.conf)
  copy_snppipeline_data.py configurationFile $scratch_out

  # Find what reads we're using. Assume all reads have been shuffled.
  READS=$(find $base_dir $base_dir/reads $base_dir/reads/shuffled -maxdepth 1 -name '*.fastq.gz' 2>/dev/null);

  # Deshuffle reads into the correct directories
  echo $READS | xargs -P $NSLOTS -n 1 sh -c '
    sample=$(basename $0 .fastq.gz);
    sampleDir=$scratch_out/samples/$sample;
    mkdir -p $sampleDir;
    echo "Deshuffling into $sampleDir";
    run_assembly_shuffleReads.pl $0 -d -gz 1>$sampleDir/1.fastq.gz 2>$sampleDir/2.fastq.gz
    if [ $? -gt 0 ]; then echo "ERROR deshuffling $i"; exit 1; fi;
  ';

  # Run snp-pipeline on the scratch drive without qsub inception.
  run_snp_pipeline.sh -c $scratch_out/snppipeline.conf -s $scratch_out/samples -m copy -o $scratch_out $REF
  if [ $? -gt 0 ]; then echo "ERROR with run_snp_pipeline.sh"; exit 1; fi;
  rm -rvf $scratch_out/samples/*

  mv -v $scratch_out $base_dir/snp-pipeline
END_OF_SCRIPT

