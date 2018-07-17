## copy this folder where you have your ROOT/BAM file

SAMPLEID=$1		## individual's ID
NAME=genome		## e.g. hiseq_pcrfree_hc_130506
FASTA=$2			## path to FASTAs
BINSIZE=100		## binsize e.g. 100 for high coverage (trio data), 1000 for low coverage (1KG data)
ROOT=$3			  ## path of (tree) ROOT file from CNVnator
SNP=$4		  	## path to SNV file
NATOR=$5      ## path to cnvnator directory

## after running CNVnator on the BAM and produces a ROOT file, this generates the histogram from the ROOT file using CNVnator and the FASTAs
ln -s $ROOT
$NATOR/cnvnator -root $ROOT -outroot his.$SAMPLEID.$NAME.root -his $BINSIZE  -d $FASTA
$NATOR/cnvnator -root his.$SAMPLEID.$NAME.root -stat $BINSIZE
$NATOR/cnvnator -root his.$SAMPLEID.$NAME.root -eval $BINSIZE > binSize$BINSIZE.log


## prepare addRD file
rd=$(grep "Average RD per bin (1-22) is" binSize$BINSIZE.log | sed 's/Average RD per bin (1-22) is //g'  | awk '{printf("%d\n"),$SAMPLEID + 0.5}')
./print_addRDcpp.sh $ROOT $rd./$BINSIZE
make addRD

## run addRD
ln -s $SNP
./addRD $SNP >& rd.$SAMPLEID.$NAME.cnvnator.log
