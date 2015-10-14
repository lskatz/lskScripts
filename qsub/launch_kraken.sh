#!/bin/bash
#$ -pe smp 1-16
#$ -S /bin/bash
#$ -cwd -V
#$ -o kraken.log -j y
#$ -N kraken

# The aspen module for kraken is broken, and so I just have
# to assume Kraken is in the path.

#source /etc/profile.d/modules.sh
#module load kraken/0.10.4

function logmsg () {
  script=$(basename $0);
  echo -e "$script: $@" >&2;
}

function run () {
  script=$(basename $0);
  logmsg "Running $@"
  eval $@
  if [ $? -gt 0 ]; then logmsg "ERROR with previous command"; fi;
}

NSLOTS=${NSLOTS-1}
KRAKEN_DEFAULT_DB=${KRAKEN_DEFAULT_DB-/scicomp/reference/kraken/0.10.4/mini-20140330}

if [ $# -eq 0 ]; then
  logmsg "Usage: $0 out.html reads_1.fastq.gz reads_2.fastq.gz [more_1.fastq.gz more_2.fastq.gz ...]";
  logmsg "NOTE: KRAKEN_DEFAULT_DB is currently set to $KRAKEN_DEFAULT_DB"
  exit 1;
fi;

# Output HTML is the first arg.
# The rest of the args are reads.
HTML=$1; shift;
READS=$@;


# Where are my executables?
KRAKENDIR=$(dirname $(which kraken));
KRONADIR=$(dirname $(which ktImportText));

# Set up the work space
TEMPDIR=$(mktemp --directory --suffix=$(basename $0));
KRAKENOUT="$TEMPDIR/kraken.out"
KRAKENTAXONOMY="$TEMPDIR/kraken.taxonomy";

logmsg "tempdir is $TEMPDIR\n  kraken dir is $KRAKENDIR\n  krona dir is $KRONADIR";

# Run the Kraken to Krona pipeline but don't exit, 
# so that TEMPDIR can be deleted each time.

run $KRAKENDIR/kraken --fastq-input --paired --db=$KRAKEN_DEFAULT_DB --preload --gzip-compressed --quick --threads $NSLOTS --output $KRAKENOUT $READS

# kraken-translate --db $KRAKEN_DEFAULT_DB --mpa-format $KRAKENOUT > translate
run kraken-translate --db $KRAKEN_DEFAULT_DB $KRAKENOUT | cut -f 2- | sort | uniq -c |\
  perl -lane '
              s/^ +//;   # remove leading spaces
              s/ +/\t/;  # change first set of spaces from uniq -c to a tab
              s/;/\t/g;  # change the semicolon-delimited taxonomy to tab-delimited
              print;
             ' |\
  sort -k1,1nr > $KRAKENTAXONOMY

# $KRAKENDIR/kraken-report --db $KRAKEN_DEFAULT_DB --show-zeros $KRAKENOUT> $KRAKENTAXONOMY.withzeros 2>/dev/null &
#run $KRAKENDIR/kraken-report --db $KRAKEN_DEFAULT_DB $KRAKENOUT> $KRAKENTAXONOMY
#logmsg "There are " $(wc -l $KRAKENTAXONOMY) " lines in the taxonomy file"

run $KRONADIR/ktImportText -o $HTML $KRAKENTAXONOMY

rm -rvf $TEMPDIR

