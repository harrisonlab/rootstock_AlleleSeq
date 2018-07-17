#!/bin/bash

#Assemble contigs using Bowtie
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

R1=$1
REF=$2
DEST=$3
OUTNAME=$4
WORK_DIR=$TMPDIR

echo  "Running Bowtie 2S with the following in= REF IS '$REF' READ 1 '$R1' READ 2 ' $R2 ' $DEST "
cd $WORK_DIR

#bowtie2 -p 4 -x $REF -U $R1  -S $WORK_DIR/$OUTNAME
bowtie2 -p 4 --no-unal --un-conc $OUTNAME -x $REF -U $R1 -S $OUTNAME


cp * $DEST/.


