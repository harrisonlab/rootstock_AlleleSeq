#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

REF=$1
BAMS=$2
DEST=$3
OUTNAME=$4
REGIONS=$5

count=0
while read REGION
do
	let count++
	qsub $SCRIPT_DIR/submit_pileup2.sh $REF $BAMS $DEST $OUTNAME.$count $REGION
done < $REGIONS 