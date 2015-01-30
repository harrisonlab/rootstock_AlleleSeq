#!/bin/bash
SCRIPT_DIR=$(readlink -f ${0%/*})
FWD
REV
TRIMDIR
	qsub $SCRIPT_DIR/submit_fastq-mcf.sh $SCRIPT_DIR $FWD $REV $TRIM_DIR

done 
