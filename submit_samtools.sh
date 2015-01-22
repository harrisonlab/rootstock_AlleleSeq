#!/bin/bash

#Assemble contigs using Bowtie
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

SAM=$1
BAM=$2
BAMOUT=$3
SORTEDOUT=$4
WORK_DIR=$TMPDIR

echo  "Running Samtools with $SAM to get $OUT and output to $BAM\n"

samtools view -S -b $SAM >$WORK_DIR/$BAMOUT
samtools sort $WORK_DIR/$BAMOUT $WORK_DIR/$SORTEDOUT


echo "Completed samtools\n"
$WORK_DIR
cp -r $WORK_DIR/* $BAM/.
echo "files copied"
