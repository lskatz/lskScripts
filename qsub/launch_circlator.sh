#!/bin/bash
#$ -e assembly.circlator.err
#$ -o assembly.circlator.out
#$ -N circ.assembly
#$ -pe smp 4-16
source /etc/profile.d/modules.sh
module load Python/3.4
module load circlator/1.2.1
if [ -d circ_dir ]; then rm -r circ_dir; fi 
circlator all --threads $NSLOTS assembly.fasta reads.fasta circ_dir 
