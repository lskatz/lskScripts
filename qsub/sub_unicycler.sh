#!/bin/bash
#$ -e err.err
#$ -o out.out
#$ -N racon
#$ -pe smp 8,16
#$ -q short.q,all.q
#$ -cwd

# Author: Dhwani Batra

source /etc/profile.d/modules.sh
module purge
unset PYTHONPATH
module load Unicycler/0.4.4



usage(){

        echo -e "\nUSAGE: $(basename $0)\n"\
                "   -1  Full path to Illumina Read 1 \n"\
                "   -2  Full path to Illumina Read 2 \n"\
		"   -r  Full path to Pacbio Reference \n"\
        exit 1

}

SHORT1=""
SHORT2=""
REFERENCE=""

while getopts 1:2:r: opt; do
        case $opt in
        1)
                SHORT1=$OPTARG
                ;;
        2)
                SHORT2=$OPTARG
                ;;
	r)		
                REFERENCE=$OPTARG
		;;

        esac
done

if [ -z "$SHORT1" ]; then
        echo -e "\nERROR: -1 is a required parameter."
        usage
        exit 1
fi

if [ -z "$SHORT2" ]; then
        echo -e "\nERROR: -2 is a required parameter."
        usage
        exit 1
fi

if [ -z "$REFERENCE" ]; then
		echo $REFERENCE
        echo -e "\nERROR: -r is a required parameter."
        usage
        exit 1
fi

if [ -z "$NSLOTS" ]; then
        NSLOTS=8
fi

LNAME=$(basename $REFERENCE)
LABEL=${LNAME%%.fsa}
mkdir $LABEL
cd $LABEL


unicycler_polish -1 ${SHORT1} -2 ${SHORT2} -a ${REFERENCE} --pilon /apps/x86_64/pilon/1.22/lib/pilon/pilon-1.22.jar --ale /apps/x86_64/ALE/ALE-20130717/src/ALE  --threads $NSLOTS
cd ..
