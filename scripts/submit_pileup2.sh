#!/bin/bash

#Pileup BAM files with samtools
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

REF=$1
#BAMS=$(awk -vRS="\n" -vORS=" " '1' $2)
BAMS=$2
DEST=$3
OUTNAME=$4
REGION=$5
WORK_DIR=$TMPDIR

echo "piling up: $BAMS"
echo "ref genome: $REF" 
echo "region: $REGION" 
echo "output file: $OUTNAME"

#cd $WORK_DIR

samtools mpileup -o $OUTNAME -B -u -t DP,DPR,DP4 -r $REGION -f $REF -b $BAMS

#cp $OUTNAME $DEST/. 