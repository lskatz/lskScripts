#!/bin/bash
# Author: Lee Katz <lkatz@cdc.gov>

#$ -S /bin/bash
#$ -pe smp 4-24
#$ -cwd -V
#$ -o chewbbaca.log
#$ -j y
#$ -N Chewbbaca

set -eu

scriptname=$(basename $0);
dirname=$(dirname $0);
if [ "$#" -lt 3 ]; then
  echo "Runs cgMLST on a genome assembly"
  echo "Usage: $scriptname outdir/ cgMLST-database/ indir"
  echo "  where indir is a directory of assembly fasta files"
  exit 0;
fi;

outdir=$1
DB=$2
indir=$3

# Some defaults
NSLOTS=${NSLOTS:=24}
TMPDIR=${TMPDIR:=/scratch}

function logmsg() {
  script=$(basename $0)
  HR="-------------"
  echo $HR >&2
  echo "$script: $@" >&2
  echo $HR >&2
}

if [ ! -d "$DB" ]; then
  logmsg "ERROR: could not find a cgMLST database folder at $DB"
  exit 1
fi

# Count how many fasta files are in the directory and if
# they aren't there then make an error
num_schema_loci=$(\ls -f1 "$DB" | grep -c 'fasta$' || true)
if [ "$num_schema_loci" -lt 1 ]; then
  logmsg "ERROR: no fasta files found in $DB"
  exit 1
fi

logmsg "Found $num_schema_loci loci in the schema at $DB"

tempdir=$(mktemp --tmpdir=$TMPDIR --directory chewbbaca.XXXXXX)
trap ' { cd; rm -rf $tempdir; } ' EXIT

# Runs chewbbaca in a container
# Arguments:
#   database path
#   file of fasta filenames or fasta file
function chewbbaca() {
  DB=$1
  indir=$2
  # TMPDIR is a global for basically /scratch
  # $tempdir is a global for a directory under /scratch
  # NSLOTS

  ls -dl $DB
  ls -dl $indir
  ls -dl $TMPDIR
  ls -dl $tempdir
  echo "NSLOTS: $NSLOTS"

  rm -rfv $DB/temp
  set -x
  singularity exec -B $TMPDIR:$TMPDIR -B $PWD:$PWD -B $indir:/input -B $DB:/schema -B $tempdir:/out bin/chewbbaca.cif chewBBACA.py AlleleCall -i /input --schema-directory /schema -o /out --cpu $NSLOTS 
  set +x
  mv -v $tempdir/results* $outdir/ || true

  logmsg $tempdir
  ls -lhA $tempdir
}

mkdir -pv "$tempdir/input"
cp -vL $indir/*.f*a $tempdir/input/ || true
#cp -vrL $indir "$tempdir/input" || true

mkdir -p $outdir
chewbbaca $DB "$tempdir/input"

