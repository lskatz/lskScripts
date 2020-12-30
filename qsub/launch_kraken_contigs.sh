#!/bin/bash
#$ -pe smp 1-16
#$ -S /bin/bash
#$ -cwd -V
#$ -o kraken.log -j y
#$ -N kraken

# The aspen module for kraken is broken, and so I just have
# to assume Kraken is in the path.

source /etc/profile.d/modules.sh
#module load kraken/0.10.4
module load kraken/1.0.0
module load krona/2.5
export PATH=$PATH:~/src/lskScripts/scripts

function logmsg () {
  script=$(basename $0);
  echo -e "$script: $@" >&2;
}

function run () {
  script=$(basename $0);
  logmsg "Running $@"
  eval $@
  if [ $? -gt 0 ]; then logmsg "ERROR with previous command"; exit 1; fi;
}

NSLOTS=${NSLOTS-8}
KRAKEN_DEFAULT_DB=${KRAKEN_DEFAULT_DB-/scicomp/reference/kraken/0.10.4/mini-20140330}

if [ $# -eq 0 ]; then
  logmsg "Usage: $0 outdir in.fasta";
  logmsg "NOTE: KRAKEN_DEFAULT_DB is currently set to $KRAKEN_DEFAULT_DB"
  exit 1;
fi;

OUTDIR=$1
ASM=$2

# Where are my executables?
KRAKENDIR=$(dirname $(which kraken));
KRONADIR=$(dirname $(which ktImportText));

# Set up the work space
if [ -d "$OUTDIR" ]; then 
  echo "$OUTDIR already exists"; 
  exit 1; 
fi;
mkdir $OUTDIR
if [ $? -gt 0 ]; then
  echo "ERROR: I could not create $OUTDIR";
  exit 1
fi
KRAKENOUT="$OUTDIR/kraken.out"
KRAKENTAXONOMY="$OUTDIR/kraken.taxonomy";
HTML="$OUTDIR/kraken.html"

logmsg "Outdir is $OUTDIR\n  kraken dir is $KRAKENDIR\n  krona dir is $KRONADIR";

run $KRAKENDIR/kraken --fasta-input --db=$KRAKEN_DEFAULT_DB --preload --threads $NSLOTS --output $KRAKENOUT $ASM

translate-kraken-contigs.pl $KRAKENOUT | sort -k1,1nr > $KRAKENTAXONOMY

run $KRONADIR/ktImportText -o $HTML $KRAKENTAXONOMY



