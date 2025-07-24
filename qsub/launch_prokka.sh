#!/bin/bash
#$ -pe smp 8-16
#$ -S /bin/bash
#$ -cwd -V
#$ -o prokka.log -j y
#$ -N prokka

contigs=$1
genome=$2

#source /etc/profile.d/modules.sh
conda activate prokka

genus=${3-genus}
species=${4-species}
NSLOTS=${NSLOTS-1}

if [ "$genome" == "" ]; then
  script=$(basename $0);
  echo "Usage: $script contigs.fasta genomename [genus species]"
  exit 1;
fi
#module load prokka/1.13.3
#module load rnammer/1.2

command="prokka --prefix $genome --compliant --locustag $genome --genus $genus --species $species --strain $genome --force --cpus $NSLOTS $contigs"
eval $command
