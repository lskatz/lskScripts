#!/bin/sh
# Author: Lee Katz <gzu2@cdc.gov>

#$ -cwd -V 
#$ -S /bin/sh
#$ -q all.q -pe smp 1
#$ -N downloadSrrViaRGN
#$ -o download.log -j y

## CONFIGURATION
RGN_ASCP_PATH="/opt/aspera/bin/ascp"
ASCP_XOPTS="-v -QT -l640M -i /opt/aspera/etc/asperaweb_id_dsa.putty"
RGN_URI="gzu2-rgntds@rgntds.cdc.gov"
## END CONFIGURATION
#####################

SRR=$1
NAME=$2
OUTDIR=$3
THISSCRIPT=$(basename $0);

# Usage statement
if [ "$NAME" == "" ]; then
  echo "Downloads a genome using the RGN, then sends it back to you, then decompresses it into split reads"
  echo "Usage: $THISSCRIPT SRR0123456 nameOfGenome outdir"
  echo "Example: $THISSCRIPT SRR1041486 FL_FLDACS-00090 ."
  exit 1;
fi


echo `date +'%H:%M:%S'`" Transferring the file to the remote computer from NCBI"
THREE=${SRR:0:3}
SIX=${SRR:0:6}
ssh $RGN_URI "mkdir -p ~/tmp; $RGN_ASCP_PATH $ASCP_XOPTS anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$THREE/$SIX/$SRR/$SRR.sra ~/tmp/$SRR.sra"
if [ $? -gt 0 ]; then echo "ERROR with ascp on the remote computer!"; exit 1; fi;

echo `date +'%H:%M:%S'`" Transferring the file back to this computer"
rsync --progress -a $RGN_URI:~/tmp/$SRR.sra $OUTDIR/$NAME.sra
if [ $? -gt 0 ]; then echo "ERROR with transferring the file to here using rsync"; exit 1; fi;

echo `date +'%H:%M:%S'` "Decompressing the file into fastq.gz - this might take a while";
fastq-dump -v --defline-seq '@$ac_$sn[_$rn]/$ri' --defline-qual '+' --split-files -O . --gzip $OUTDIR/$NAME.sra
if [ $? -gt 0 ]; then echo "ERROR with fastq-dump"; exit 1; fi;


echo `date +'%H:%M:%S'`" Finished. Files will be found in $OUTDIR";

exit 0;

