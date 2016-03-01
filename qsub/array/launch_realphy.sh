#!/bin/bash -l

# Runs any set of reads through KSNP in a cluster-friendly way.
# Each reads directory will be a distinct project.
# Author: Lee Katz
# Usage: bash $0 readsdir1 [readsdir2 ... readsdirX]

if [ "$2" == "" ]; then
  echo "Usage: $0 ref.fasta cov1 [cov2...]"
  exit 1;
fi

TMP=$(mktemp --tmpdir='.' --directory qsubRealPhy.XXXXXXXX)
echo "tmp dir is $TMP "

export REF=$1; shift;
CTRL_FILE="$TMP/array.txt"
echo "$@" | tr ' ' '\n' | grep . > $CTRL_FILE

mkdir -p $TMP/log
#qsub -q all.q -N RealPhy -o $TMP/log -j y -pe smp 12 -hard -l exclusive=true -V -S /bin/bash -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
qsub -q all.q -N RealPhy -o $TMP/log -j y -pe smp 2 -V -S /bin/bash -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "REF=$REF" -v "CTRL_FILE=$CTRL_FILE" <<- "END_OF_SCRIPT"
  #!/bin/bash

  module load RealPhy/v112
  module load bowtie2/2.2.4
  module load phylip/3.69
  module load phyml/3.0
  module load tree-puzzle/5.3.rc16

  base_dir=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE)
  echo "Working on $base_dir, host " $(hostname)
  mkdir -p /scratch/$USER
  scratch_out=$(mktemp --tmpdir="/scratch/$USER" --directory realphy.XXXXXX)
  export scratch_out
  echo "Temporary directory is $scratch_out"
  if [ ! -e "$scratch_out" ]; then echo "ERROR: could make temporary directory $scratch_out"; exit 1; fi;
  if [ -e $base_dir/RealPhy ]; then echo "Found $base_dir/RealPhy. Skipping."; exit 0; fi

  # Find what reads we're using. Do not symlink
  # because realphy has a stupid problem with symlinks.
  # Instead, copy.
  READS=$(ls $base_dir/reads/*.fastq.gz);
  mkdir -p $scratch_out/samples
  cp -v $READS $scratch_out/samples/
  cp -v $REF $scratch_out/samples/reference.fasta

  mkdir -p $scratch_out/out

  # Make a config file.  RealPhy is intolerant of
  # beginning indentation for each line and of a
  # beginning empty line.
  echo -e "BOWTIE2\t/apps/x86_64/bowtie2/bowtie2-2.2.4/bowtie2
BOWTIE2BUILDER\t/apps/x86_64/bowtie2/bowtie2-2.2.4/bowtie2-build-l
TREEPUZZLE\t/apps/x86_64/tree-puzzle/5.3.rc16/bin/puzzle
RAXML\t/scicomp/groups/OID/NCEZID/DFWED/EDLB/share/bin/lyve-SET-v1.1.4/lib/standard-RAxML-8.1.16/raxmlHPC-PTHREADS
Rscript\t/apps/x86_64/R/3.2.3/bin/Rscript
MaxPars\t/apps/x86_64/phylip/phylip-3.69/exe/dnapars
PhyML\t/apps/x86_64/phyml/PhyML_3.0/phyml
" > $scratch_out/out/config.txt

  # Sanity check: look at the files and their permissions in
  # the config file.
  cut -f 2 $scratch_out/out/config.txt | xargs -I {} ls -lh {}

  echo "Config file contents:"
  cat $scratch_out/out/config.txt

  REALPHY_v112 $scratch_out/samples $scratch_out/out -readLength 250 -ref reference
  if [ $? -gt 0 ]; then exit 1; fi;

  ls -lh $scratch_out/reference/alignOut_NoGenes > $scratch_out/reference/alignOut_NoGenes.txt
  rm -rvf $scratch_out/reference/alignOut_NoGenes
  mv -v $scratch_out/out $base_dir/RealPhy
  rm -rvf $scratch_out
END_OF_SCRIPT

