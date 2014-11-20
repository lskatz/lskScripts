#!/bin/bash
#$ -pe smp 8-16
#$ -S /bin/bash
#$ -cwd -V
#$ -o prokka.log -j y
#$ -N prokka

contigs=$1
genome=$2

genus=${3-genus}
species=${4-species}
NSLOTS=${NSLOTS-1}
command="prokka --prefix $genome --addgenes --locustag $genome --genus $genus --species $species --strain $genome --force --cpus $NSLOTS $contigs"
$command
