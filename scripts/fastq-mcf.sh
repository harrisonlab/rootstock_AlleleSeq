#!/bin/bash
SCRIPT_DIR=$(readlink -f ${0%/*})
FWD=$1
REV=$2
TRIMDIR=$3
	qsub $SCRIPT_DIR/submit_fastq-mcf.sh $SCRIPT_DIR $FWD $REV $TRIM_DIR

