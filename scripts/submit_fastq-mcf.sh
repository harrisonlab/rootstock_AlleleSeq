#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=4G

WORK_DIR=$TMPDIR
SCRIPT_DIR=$1
FORWARD_READ=$2
REVERSE_READ=$3
TRIMDIR=$4

#FASTQ MCF COMMAND
cd $WORK_DIR
fastq-mcf $SCRIPT_DIR/illumina_full_adapters.fa $FORWARD_READ $REVERSE_READ -o $FORWARD_READ.trim -o $REVERSE_READ.trim -C 10000000 -u -k 20 -t 0.01 -p 20

cp * $TRIMDIR/.
