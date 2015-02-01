#!/bin/bash

#Assemble contigs using Bowtie
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

R1=$1
R2=$2
REF=$3
DEST=$4
OUTNAME=$5
X=$6
I=$7
WORK_DIR=$TMPDIR

echo  "Running Bowtie 2S with the following in= REF IS '$REF' READ 1 '$R1' READ 2 ' $R2 ' DEST is '$DEST' OUTNAME IS '$OUTNAME' MIN is '$X' MAX is '$I'"

bowtie2 -p 4 --no-unal --un-conc $WORK_DIR/$OUTNAME -X $X -I $I -x $REF -1 $R1 -2 $R2 -S $WORK_DIR/$OUTNAME

cd $WORK_DIR
cp * $DEST/.

