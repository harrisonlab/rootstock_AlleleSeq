#!/bin/bash
#AlleleSeq pipeline
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

BASE=$1
RNA_1=$2
RNA_2=$3
PI=$4
MI=$5
SNPS=$6
CNVS=$7
BNDS=hits.bed
MAPS=$8
OUTFILE=$9
DEST=$BASE

FDR_SIMS=5
FDR_CUTOFF=0.1

PL=/home/deakig/projects/apple_rootstock/scripts/AlleleSeq
TCO=$RNA_1.out

export PYTHONPATH=~/usr/local/lib/python2.7/site-packages/

echo "Running..." 
	
bash -c  "python $PL/MergeBowtie.py \
           <($PL/filter_input.sh $PL $RNA_1 $RNA_2 a1; bowtie2 --no-hd --no-unal -f -x $PI -1 a1_1 -2 a1_2|grep -v XS |cut -f 1,2,3,4,10) \
           <($PL/filter_input.sh $PL $RNA_1 $RNA_2 a2; bowtie2 --no-hd --no-unal -f -x $MI -1 a2_1 -2 a2_2|grep -v XS |cut -f 1,2,3,4,10) \
           $MAPS | python $PL/SnpCounts.py $SNPS - $MAPS $TCO"

python $PL/CombineSnpCounts.py 5 $SNPS $BNDS $CNVS $RNA_1.counts.txt $RNA_1.counts.log $TCO
python $PL/FalsePos.py $RNA_1.counts.txt $FDR_SIMS $FDR_CUTOFF > $RNA_1.FDR.txt
awk -f $PL/finalFilter.awk thresh=$(awk 'END {print $6}' $RNA_1.FDR.txt) < $RNA_1.counts.txt > $OUTFILE

cp $OUTFILE $DEST/.
