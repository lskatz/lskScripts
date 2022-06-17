#!/bin/bash -l

# Runs NCBI's wgmlst in an array for MLST calling

OUTDIR=$1; shift;
DB=$1; shift
ASM=$@

set -e
set -u

if [ "$ASM" == "" ]; then
  echo "Run NCBI's wgmlst for wgMLST allele calling"
  echo "Usage: $0 outdir something.scheme *.fasta"
  exit 1;
fi

if [ ! -d "$OUTDIR" ]; then
  mkdir "$OUTDIR"
fi;
if [ ! -e "$DB" ]; then
  echo "ERROR: not found: $DB"
  exit 2;
fi

tmpdir=$(mktemp --tmpdir='.' --directory wgmlst.XXXXXXXX)
mkdir -pv $tmpdir/log
#trap " rm -rf $tmpdir " EXIT
echo "tmp dir is $tmpdir "

# CTRL file will have per line:
#   filename  coverageLevel
CTRL_FILE="$tmpdir/array.txt"
# Put the reads one at a time into the CTRL_FILE but use paste to keep paired ends together
echo "$ASM" | tr ' ' '\n' > $CTRL_FILE

head $CTRL_FILE

module purge

# Check executables
which wgmlst

qsub -N ncbi_wgmlst -o $tmpdir/log -j y -pe smp 1 -V -cwd -t 1-$(cat $CTRL_FILE | wc -l) \
  -v "tmpdir=$tmpdir" -v "DB=$DB" -v "CTRL_FILE=$CTRL_FILE" -v "OUTDIR=$OUTDIR" <<- "END_OF_SCRIPT"
  #!/bin/bash -l
  set -e
  set -u

  # Set up filenames
  fasta=$(sed -n ${SGE_TASK_ID}p $CTRL_FILE);
  name=$(basename $fasta .fasta)

  sampledir=$(mktemp --tmpdir=$tmpdir --directory $name.wgmlst.XXXXXX);
  trap "rm -rf $sampledir" EXIT
  finalout="$OUTDIR/$name"

  if [ -e "$finalout/.done" ]; then
    echo "Found $finalout/.done; not repeating."
    exit 0;
  fi

  mappings=$sampledir/mappings
  alleles=$sampledir/alleles
  stdout=$sampledir/wgmlst.out
  log=$sampledir/wgmlst.log

  (
    date # mark how long this takes
    set -x
    wgmlst --genome $fasta --alleles $DB --cores $NSLOTS --kmer 15 --output_mappings $mappings --output_loci $alleles 
    set +x
    date
  ) 1>$stdout 2>$log

  ls -lh $sampledir
  mv -v $sampledir $finalout
  touch $finalout/.done

END_OF_SCRIPT

