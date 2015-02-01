#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

R1=$1
REF=$3
WORK_DIR=$3
OUT=$4
echo "Assemble $REF with Bowtie $R1 \n $R2 \n $REF \n $NAME \n"


#MPI WITH QSUB 

echo "R1 is $R1"

echo "REF is $2"
echo "WORK DIR is $3"
echo "OUT is $4"
qsub $SCRIPT_DIR/submit_bowtie_se.sh $R1 $REF $WORK_DIR $OUT 

