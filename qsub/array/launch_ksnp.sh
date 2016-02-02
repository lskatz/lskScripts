#!/bin/bash -l

# Runs any set of reads through KSNP in a cluster-friendly way.
# Each reads directory will be a distinct project.
# Author: Lee Katz
# Usage: bash $0 readsdir1 [readsdir2 ... readsdirX]

TMP=$(mktemp --tmpdir='.' --directory qsubKsnp3.XXXXXXXX)
echo "tmp dir is $TMP "

CTRL_FILE="$TMP/array.txt"
echo "$@" | tr ' ' '\n' | grep . > $CTRL_FILE

mkdir -p $TMP/log
qsub -q all.q -N KSNP3 -o $TMP/log -j y -pe smp 3-4 -V -S /bin/bash -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash

  module load kSNP/3.0.0
  module load fastx-toolkit/0.0.13

  base_dir=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  echo "Working on $base_dir, host " $(hostname)
  mkdir -p /scratch/$USER
  scratch_out=$(mktemp --tmpdir="/scratch/$USER" --directory ksnp3.XXXXXX)
  export scratch_out
  echo "Temporary directory is $scratch_out"
  if [ ! -e "$scratch_out" ]; then echo "ERROR: could make temporary directory $scratch_out"; exit 1; fi;

  # Find what reads we're using
  READS=$(find $base_dir $base_dir/reads $base_dir/reads/shuffled -maxdepth 1 -name '*.fastq.gz');
  echo -e "Found reads\n" $READS

  # Convert to fasta
  echo $READS | xargs -P $NSLOTS -n 1 sh -c '
    sample=$(basename $0 .fastq.gz);
    sampleDir=$scratch_out/samples/$sample;
    mkdir -pv $sampleDir;

    zcat $0 | fastq_to_fasta -Q33 > $sampleDir/$sample.fasta
    if [ $? -gt 0 ]; then echo "ERROR converting $0 to fasta"; exit 1; fi;
    merge_fasta_reads3 $sampleDir/$sample.fasta > $sampleDir/merged.fasta
    if [ $? -gt 0 ]; then echo "ERROR merging $0"; exit 1; fi;

    # cleanup
    mv $sampleDir/merged.fasta $scratch_out/samples/$sample.fasta
    rm -rvf $sampleDir
  ';

  # Switch to the tmp dir
  pushd $scratch_out

    # Find optimal kmer value
    MakeKSNP3infile samples in.txt A
    cat samples/*.fasta > tmp.fasta
    KMERLENGTH=$(Kchooser tmp.fasta | grep "The optimum value of K is" | grep -o "[0-9]\+")
    if [ $? -gt 0 ]; then echo "ERROR with Kchooser"; exit 1; fi;
    echo "The optimal kmer length is $KMERLENGTH";
    rm -f tmp.fasta

    kSNP3 -k $KMERLENGTH -in in.txt -c 4 -core -ML -min_frac 0.75 -CPU $NSLOTS -NJ -ML -vcf -outdir out
    if [ $? -gt 0 ]; then echo "ERROR with kSNP3"; exit 1; fi;

  # pop out of the tmp dir
  popd

  mv -v $scratch_out/out $base_dir/kSNP3
  rm -rvf $scratch_out
END_OF_SCRIPT

