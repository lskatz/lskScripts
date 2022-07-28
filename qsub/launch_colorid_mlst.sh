#!/bin/bash
# Runs ColorID MLST
# Author: Lee Katz <lkatz@cdc.gov>

#$ -S /bin/bash
#$ -pe smp 1
#$ -cwd -V
#$ -o EToKi.log
#$ -j y
#$ -N EToKi

outdir=$1
dbdir=$2
asm=$3

NSLOTS=${NSLOTS:=1}

source /etc/profile.d/modules.sh
scriptname=$(basename $0);


if [ "$asm" == "" ]; then
  echo "Usage: $scriptname outdir dbdir asm.fasta"
  echo "  where dbdir is a directory of locus fasta files like chewbbaca"
  exit 0;
fi;

set -e
set -u

module purge

colorid=$(which colorid_bv || which colorid)
process_MLST=$(which process_MLST.py)

tmpdir=$(mktemp --directory --suffix=$(basename $0));
trap ' { rm -rf $tmpdir; } ' EXIT

mkdir -pv $outdir || echo "WARNING: outdir already exists: $outdir"

bxi="$tmpdir/asm.bxi"

# Estimate assembly size by file size and then multiply by 10x
asm_size=$(cat $asm | wc -c)
bxi_size=$(( $asm_size * 10 ))

# Some colorid vars
fofn="$tmpdir/asm.fofn"
alleles="$tmpdir/alleles.tsv"

# Get a file of filenames
echo -e "$(basename $asm .fasta)\t$asm" > $fofn

$colorid build -b $bxi -s $bxi_size -n 2 -k 39 -t $NSLOTS -r $fofn

$colorid search -ms -b $bxi.bxi -q $dbdir/*.fasta > $alleles
sed -i.bak '/\*/d' $alleles
$process_MLST $alleles $tmpdir/mlst

cp -v $alleles $tmpdir/mlst* $outdir/

