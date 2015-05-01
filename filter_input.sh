# This script should take whatever input is provided and stream a fasta format suitable for bowtie2

PL=$1
INFILE1=$2
INFILE2=$3
q=$4

if [ -p ${q}_1 ]; then
	rm ${q}_1
	mkfifo ${q}_1
else
	mkfifo ${q}_1 
fi

if [ -p ${q}_2 ]; then
	rm ${q}_2
	mkfifo ${q}_2
else
	mkfifo ${q}_2 
fi


cat ${INFILE1}|python ${PL}/fastq2result.py - -|python ${PL}/ConvertTags.py >${q}_1 & 
cat ${INFILE2}|python ${PL}/fastq2result.py - -|python ${PL}/ConvertTags.py >${q}_2 &
