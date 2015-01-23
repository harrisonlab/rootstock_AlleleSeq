#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

SAM=$1
BAMDIR=$2
OUT=$3
OUT2=$4
echo "Take $SAM and make to $OUT and then sort and output to $OUT2 \n"
echo "SAM is $SAM"
echo "BAMDIR is $BAMDIR"
echo "BAM is $OUT"
echo "SORTED is $OUT2"

qsub $SCRIPT_DIR/submit_samtools.sh $SAM $BAMDIR $OUT $OUT2

