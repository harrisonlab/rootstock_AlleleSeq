#!/bin/bash

#Samtools 
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

SAM=$1
BAMDIR=$2
BAMOUT=$3
SORTEDOUT=$4
WORK_DIR=$TMPDIR

echo  "Running Samtools with $SAM"
echo  "to get $BAMOUT "
echo  "and output to $BAMDIR\n"
echo "THEN SORTING AND OUTPUTTING TO $SORTEDOUT \n"
echo "SORTING THIS  $BAMDIR/$BAMOUT "
echo "OUTPUTTING THIS  $BAMDIR/$SORTEDOUT "

#which samtools

samtools view -S -b $SAM >$WORK_DIR/$BAMOUT
samtools sort $WORK_DIR/$BAMOUT $WORK_DIR/$SORTEDOUT

cd $WORK_DIR
cp * $BAMDIR/.

