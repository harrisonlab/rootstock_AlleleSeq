#!/bin/bash

SCRIPT_DIR=$(readlink -f ${0%/*})

R1=$1
R2=$2
REF=$3
WORK_DIR=$4
echo "Assemble $REF with Bowtie $R1 \n $R2 \n $REF \n $NAME \n"

#MPI WITH QSUB 

echo "R1 is $R1"
echo "R2 is $R2"
echo "REF is $REF"
echo "WORK DIR is $4"

qsub $SCRIPT_DIR/submit_bowtie.sh $R1 $R2 $REF $WORK_DIR

