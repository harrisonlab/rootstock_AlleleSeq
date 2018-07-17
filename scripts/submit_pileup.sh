#!/bin/bash

#Pileup BAM files with samtools
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

REF=$1
BAMS=$(awk -vRS="\n" -vORS=" " '1' $2)
DEST=$3
OUTNAME=$4
WORK_DIR=$TMPDIR

echo "piling up $BAMS with $REF output file: $OUTNAME"

cd $WORK_DIR

samtools mpileup -o $OUTNAME -uf $REF $BAMS

cp * $DEST/.
