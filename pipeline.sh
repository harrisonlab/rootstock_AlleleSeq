#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

BASE=$1
RNA_1=$2
RNA_2=$3
PI=$BASE/pat2/$4
MI=$BASE/mat2/$5
SNPS=$BASE/$6
CNVS=$BASE/$7
BNDS=hits.bed
MAPS=$BASE/chr%s_$8
OUTFILE=$9



qsub $SCRIPT_DIR/submit_pipeline.sh $BASE $RNA_1 $RNA_2 $PI $MI $SNPS $CNVS $MAPS $OUTFILE
