#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

SAM=$1
BAMDIR=$2
OUT=$3
echo "Take $SAM and make to $BAM and then sort and output to $OUT \n"

#MPI WITH QSUB

echo "SAM is $SAM"
echo "BAMDIR is $BAMDIR"
echo "BAM is $OUT"

qsub $SCRIPT_DIR/submit_samtools.sh $SAM $BAMDIR $OUT

