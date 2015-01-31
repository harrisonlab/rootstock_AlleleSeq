#!/bin/bash

#Assemble contigs using Bowtie
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

R1=$1
REF=$2
DEST=$3
WORK_DIR=$TMPDIR

echo  "Running Bowtie 2S with the following in= REF IS '$REF' READ 1 '$R1' READ 2 ' $R2 ' $DEST "

bowtie2 -p 4 -x $REF -U $R1  -S $DEST

$WORK_DIR
#cp -r $WORK_DIR/* $DEST/.
#echo "files copied"
