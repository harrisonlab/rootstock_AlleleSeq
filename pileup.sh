#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

REF=$1
BAMS=$2
DEST=$3
OUTNAME=$4

qsub $SCRIPT_DIR/submit_pileup.sh $REF $BAMS $DEST $OUTNAME 
