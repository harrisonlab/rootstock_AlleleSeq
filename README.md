# rootstock

To run this analysis I used:
fastqc
bowtie2
samtools
bcftools
beagle

$ROOTSTOCK was set to the directory containg the rootstock project (~/projects/apple_rootstock/rootstock_genetics)

Folders for each rootstock was created in the rootstock_genetics folder

Data was unzipped and concatenated into R1 and R2 files in the subdirectory conc within each rootstock folder
eg:
```shell
gunzip *
cat 863_LIB6292_LDI5172_GTGAAA_L00*_R1.fastq.gz >m27_r1.fa
```

within the conc directory of each genome folder, fastqc was run on each pair of reads in order to assess the quality
```shell
nohup fastqc m116_r1.fq m116_r2.fq &
nohup fastqc m27_r1.fq m27_r2.fq &
nohup fastqc m9_r1.fq m9_r2.fq &
nohup fastqc m13_r1.fq m13_r2.fq &
nohup fastqc mm106_r1.fq mm106_r2.fq &
nohup fastqc gala_r1.fq gala_r2.fq &
nohup fastqc o3_r1.fq o3_r2.fq &

 ```

Trimming was performed with fastq-mcf (o3 is SE)
```shell
cd $ROOTSTOCK/rootstock_genetics

for d in [mg]*;
 do
  $ROOTSTOCK/scripts/fastq-mcf.sh $ROOTSTOCK/rootstock_genetics/$d/conc/${d}_r1.fq $ROOTSTOCK/rootstock_genetics/$d/conc/${d}_r2.fq $ROOTSTOCK/rootstock_genetics/$d/conc/  
 done

$ROOTSTOCK/scripts//fastq-mcf.sh $ROOTSTOCK/rootstock_genetics/o3/conc/o3_r1.fq  $ROOTSTOCK/rootstock_genetics/o3/conc/ 
 ```
 
 Removal of phix
 ```shell
 cd $ROOTSTOCK/rootstock_genetics

for d in [mg]*;
 do
  $ROOTSTOCK/scripts/bowtie.sh $ROOTSTOCK/rootstock_genetics/$d/conc/${d}_r1.fq.trim $ROOTSTOCK/rootstock_genetics/$d/conc/${d}_r2.fq.trim /home/groups/harrisonlab/ref_genomes/phix/phix $ROOTSTOCK/rootstock_genetics/$d/conc/ phix_filtered 250 500 
 done

./bowtie_se.sh $ROOTSTOCK/rootstock_genetics/o3/conc/o3_r1.fq.trim /home/groups/harrisonlab/ref_genomes/phix/phix $ROOTSTOCK/rootstock_genetics/o3/conc/ phix_filtered 

 ```
  
##Assembly to reference
A hash was created for version 1 of the genome

```shell
cd /home/groups/harrisonlab/project_files/rootstock_genetics/
cd ref
mkdir v1
cd v1
bowtie2-build Malus_x_domestica.v1.0-primary.pseudo.fa Md
```

Made a shell script to submit automatically to grid engine
for version 1 of the genome!!

```shell
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m27/conc/phix_filtered.1 /home/groups/harrisonlab/project_files/rootstock_genetics/m27/conc/phix_filtered.2 /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m27/analysis_v1/ m27_v1 250 500
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m116/conc/phix_filtered.1  /home/groups/harrisonlab/project_files/rootstock_genetics/m116/conc/phix_filtered.2 /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m116/analysis_v1/ m116_v1 250 500
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m9/conc/phix_filtered.1  /home/groups/harrisonlab/project_files/rootstock_genetics/m9/conc/phix_filtered.2 /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m9/analysis_v1/ m9_v1 250 500
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m13/conc/phix_filtered.1  /home/groups/harrisonlab/project_files/rootstock_genetics/m13/conc/phix_filtered.2 /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m13/analysis_v1/ m13_v1 250 500
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/conc/phix_filtered.1  /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/conc/phix_filtered.2 /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/analysis_v1/ mm106_v1 250 500
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/gala/conc/phix_filtered.1  /home/groups/harrisonlab/project_files/rootstock_genetics/gala/conc/phix_filtered.2 /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/groups/harrisonlab/project_files/rootstock_genetics/gala/analysis_v1/ gala_v1 250 500

./bowtie_se.sh /home/groups/harrisonlab/project_files/rootstock_genetics/o3/conc/phix_filtered.1  /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/groups/harrisonlab/project_files/rootstock_genetics/o3/analysis_v1/ o3_v1.sam 

./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/test/conc/phix_filtered.1 /home/groups/harrisonlab/project_files/rootstock_genetics/test/conc/phix_filtered.2 /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/groups/harrisonlab/project_files/rootstock_genetics/test/analysis_v1/ test_v1 250 500

```

SAM to BAM and sort
```shell
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m27/analysis/m27_v1.sam /home/groups/harrisonlab/project_files/rootstock_genetics/m27/analysis/ m27_v1.bam m27_v1.sorted 
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m9/analysis/m9_v1.sam /home/groups/harrisonlab/project_files/rootstock_genetics/m9/analysis/ m9_v1.bam m9_v1.sorted
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m116/analysis/m116_v1.sam /home/groups/harrisonlab/project_files/rootstock_genetics/m116/analysis/ m116_v1.bam m116_v1.sorted
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m13/analysis/m13_v1.sam /home/groups/harrisonlab/project_files/rootstock_genetics/m13/analysis/ m13_v1.bam m13_v1.sorted 
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/analysis/mm106_v1.sam /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/analysis/ mm106_v1.bam mm106_v1.sorted
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/gala/analysis/gala_v1.sam /home/groups/harrisonlab/project_files/rootstock_genetics/gala/analysis/ gala_v1.bam gala_v1.sorted
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/o3/analysis/o3_v1.sam /home/groups/harrisonlab/project_files/rootstock_genetics/o3/analysis/ o3_v1.bam o3_v1.sorted

./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/test/analysis_v1/test_v1 /home/groups/harrisonlab/project_files/rootstock_genetics/test/analysis_v1 test.bam test.sorted

```

Note- to index and then view the bam file using a simple text viewer
samtools index m9.sorted.bam
samtools tview m9.sorted.bam ../../ref/v1/Malus_x_domestica.v1.0-primary.pseudo.fa

Once indexed the program qualimap can be used with the BAM files to view statistics such as coverage and insert size etc

Index the reference genome for SAMTOOLS

```shell
cd /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1
samtools faidx Malus_x_domestica.v1.0-primary.pseudo.fa 

```

Pileup into a single vcf with v1 (http://biobits.org/samtools_primer.html)

mpileup is single threaded. Multiple instances can be launched to pileup different chromosome regions. 
First step is to index all BAM files (included in bam_files txt file). Then run pileup2.sh
regions have been defined in the regions file. Finally concatenate output files (the regions file is sorted alphabetically rather than by chromosome - temp fix is to create an ordered list of output file names)
```shell
cat bam_files|xargs -I file samtools index file 
./pileup2.sh /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Malus_x_domestica.v1.0-primary.pseudo.fa bam_files /home/groups/harrisonlab/project_files/rootstock_genetics piledup.bcf regions

bcftools concat -f files >all_piledup.bcf
```
Filter output for variants 

```shell
bcftools call -Ov -v -m all_piledup.bcf > all_piledup.vcf
cat all_piledup.vcf|vcfutils.pl varFilter -d100 > flt_all_piledup.vcf
```

Define the pedigree for beagle
1) pedigree ID, 2) individual ID, 3) father’s ID, and 4) mother’s ID

m116 m116 mm106 m27
mm106 mm106 0 0
m27 m27 m13 m9
m13 m13 0 0
m9 m9 0 0
o3 o3 0 0

NB- To generate a whole consensus sequence
samtools mpileup -uf ref.fa aln.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq i 

##Phasing with Beagle
http://faculty.washington.edu/browning/beagle_utilities/utilities.html
Note SHAPE could also be used- one advantage here is it can be 'read aware'  https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#readaware

NOTE: The beagle pedigree file does not conform with the pedigree file Linkage-format and requires the pedigree ID and individual ID to be the same.

```shell

 java -jar $ROOTSTOCK/scripts/beagle.r1399.jar gt=corrected_sort_root.vcf ped=pedigree.ped out=beagle phase-its=8 overlap=1000000 window=5000000  ibd=true


#java –Xmx 12000m –jar beagle.jar gt=./beagle/flt_all_piledup.vcf ped=./beagle/pedigree.ped out=beagle chrom=1 nthreads=4
#samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf
#bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf
#java –Xmx 12000m –jar beagle.jar gt=./beagle/var.flt.vcf ped=./beagle/pedigree.ped out=beagle chrom=1 nthreads=4

```
### Phasing without beagle
Beagle does not phase as well as we would like - or the data is poor..

phasing.r can phase m27 and m116 given the list of snps as input:
```shell
# create list of snps for phasing
echo -e chr'\t'pos'\t'ref'\t'alt'\t'm9_p'\t'm9_m'\t'm27_p'\t'm27_m'\t'm116_p'\t'm116_m'\t'mm106_p'\t'mm106_m'\t'm13_p'\t'm13_m > ~/projects/apple_rootstock/allele/figure_1/snps

grep ^chr ~/projects/apple_rootstock/rootstock_genetics/beagle/corrected__snp_sort_root.vcf|awk -F"\t" '{print $1,$2,$4,$5,$11,$12,$13,$14}' OFS="\t" |perl -pe 's/:(.*?)\t/\t/g'|sed 's/:.*//g'|awk -F"/" '{$NF=$NF;print $0}' OFS="\t" >> ~/projects/apple_rootstock/allele/figure_1/snps
```

A VCF file can be produced from the output to be run through the rest of the pipeline


### Make diploid genomes from VCF

Downloaded vcf2diploid
requires that output directories are created before running command
Using local copy of Java 1.7 

##### pipeline samples
```shell

mkdir $ROOTSTOCK/allele/m27
mkdir $ROOTSTOCK/allele/m116

~/usr/bin/java -jar $ROOTSTOCK/vcf2diploid/vcf2diploid.jar -id m27 -chr ~/Data/apple/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf beagle.vcf -outDir $ROOTSTOCK/allele/m27
~/usr/bin/java -jar $ROOTSTOCK/vcf2diploid/vcf2diploid.jar -id m116 -chr ~/Data/apple/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf beagle.vcf -outDir $ROOTSTOCK/allele/m116

```
##### All other samples
```shell
mkdir $ROOTSTOCK/allele/m9
mkdir $ROOTSTOCK/allele/m13
mkdir $ROOTSTOCK/allele/m106
mkdir $ROOTSTOCK/allele/gala
mkdir $ROOTSTOCK/allele/03

~/usr/bin/java -jar $ROOTSTOCK/vcf2diploid/vcf2diploid.jar -id m9 -chr ~/Data/apple/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf beagle.vcf -outDir $ROOTSTOCK/allele/m9
~/usr/bin/java -jar $ROOTSTOCK/vcf2diploid/vcf2diploid.jar -id m13 -chr ~/Data/apple/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf beagle.vcf -outDir $ROOTSTOCK/allele/m13
~/usr/bin/java -jar $ROOTSTOCK/vcf2diploid/vcf2diploid.jar -id mm106 -chr ~/Data/apple/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf beagle.vcf -outDir $ROOTSTOCK/allele/m106
~/usr/bin/java -jar $ROOTSTOCK/vcf2diploid/vcf2diploid.jar -id gala -chr ~/Data/apple/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf beagle.vcf -outDir $ROOTSTOCK/allele/gala
~/usr/bin/java -jar $ROOTSTOCK/vcf2diploid/vcf2diploid.jar -id o3 -chr ~/Data/apple/Malus_x_domestica.v1.0-primary.pseudo.fa -vcf beagle.vcf -outDir $ROOTSTOCK/allele/o3
```

### Get variants from VCF

Run vcf2snp on beagle output (or other vcf file) to create list of snps
-c flag searches vcf file for given sample name

##### pipeline samples
```shell
$ROOTSTOCK/AlleleSeq_pipeline_v1.2/vcf2snp $ROOTSTOCKS/rootstock_genetics/beagle/beagle.vcf -c m27 > $ROOTSTOCK/allele/m27/m27.snv
$ROOTSTOCK/AlleleSeq_pipeline_v1.2/vcf2snp $ROOTSTOCKS/rootstock_genetics/beagle/beagle.vcf -c m116 > $ROOTSTOCK/allele/m116/m116.snv
```

##### all other samples
```shell
$ROOTSTOCK/AlleleSeq_pipeline_v1.2/vcf2snp $ROOTSTOCKS/rootstock_genetics/beagle/beagle.vcf -c mm106 > $ROOTSTOCK/allele/mm106/mm106.snv

$ROOTSTOCK/AlleleSeq_pipeline_v1.2/vcf2snp $ROOTSTOCKS/rootstock_genetics/beagle/beagle.vcf -c m13 > $ROOTSTOCK/allele/m13/m13.snv
$ROOTSTOCK/AlleleSeq_pipeline_v1.2/vcf2snp $ROOTSTOCKS/rootstock_genetics/beagle/beagle.vcf -c m9 > $ROOTSTOCK/allele/m9/m9.snv
```

### Make maternal and paternal indices

Maternal and paternal chromosomes will require indexing with bowtie for downstream steps.

```shell
mkdir $ROOTSTOCK/allele/m116/pat
cd $ROOTSTOCK/allele/m116/pat
ln -s ../*paternal.fa .
ref=""
for f in *.fa;do ref=$(echo $ref,$f); done
ref=$(echo ${ref#?})
bowtie2-build $(echo $ref) m116_paternal_index

mkdir $ROOTSTOCK/allele/m116/mat
cd $ROOTSTOCK/allele/m116/mat
ln -s ../*maternal.fa .
ref=""
for f in *.fa;do ref=$(echo $ref,$f); done
ref=$(echo ${ref#?})
bowtie2-build $(echo $ref) m116_maternal_index

mkdir $ROOTSTOCK/allele/m27/pat
cd $ROOTSTOCK/allele/m27/pat
ln -s ../*paternal.fa .
ref=""
for f in *.fa;do ref=$(echo $ref,$f); done
ref=$(echo ${ref#?})
bowtie2-build $(echo $ref) m27_paternal_index

mkdir $ROOTSTOCK/allele/m27/mat
cd $ROOTSTOCK/allele/m27/mat
ln -s ../*maternal.fa .
ref=""
for f in *.fa;do ref=$(echo $ref,$f); done
ref=$(echo ${ref#?})
bowtie2-build $(echo $ref) m27_maternal_index
```

### Find CNVs
Using  CNVnator. This requires ROOT (https://root.cern.ch).

##### ROOT install
V 6.x won't intall on current cluster due to incorrect gcc version.
Installed with version5. 

```shell
./configure linuxx8664gcc --prefix=/home/deakig/usr/local --etcdir=/home/usr/local/etc --with-python-libdir=/usr/lib --with-python-incdir=/usr/include --enable-builtin-pcre --disable-x11 
make
make install

```

UPDATE - root should be installed using a local copy of pcre NOT the root (boom boom) installed version. It is not on the worker nodes and otherwise downstream processes (CNVnator) will not work on these nodes. 
UPDATE2 - root reconfigured. For some reason the install ignores the prefix flag for the etc directory and tries to stick the files into /usr/local/etc, anyway the --etcdir flag can set this explicitly.

##### ROOTSYS (run each time)

The below sets ROOTSYS for each session (required for pipeline to work correctly)

```shell
. ~/usr/local/bin/thisroot.sh
```

##### CNVnator
Ran CNVnator with below:
 
```shell
$ROOTSTOCK/CNVnator_v0.3/src/cnvnator -root $ROOTSTOCK/allele/m116/m116.root -tree $ROOTSTOCK/rootstock_genetics/m116/analysis_v1/m116_v1.sorted.bam
$ROOTSTOCK/CNVnator_v0.3/src/cnvnator -root $ROOTSTOCK/allele/m27/m27.root -tree $ROOTSTOCK/rootstock_genetics/m27/analysis_v1/m27_v1.sorted.bam
```

##### Convert to AlleleSeq format

Output is in ROOT format which needs to be converted to an AlleleSeq CNV format*. alleleSeq_cnvScript can do this.
The wrapper script was edited to reflect our environment and the makefile needed editing to reflect the locaction of root (root puts its files in a root folder under lib/include). 

This pipeline makes a build specific to each root file, therefore the following files were copied from $ROOTSTOCK/alleleSeq_cnvScript to $ROOTSTOCK/allele/m27/m27_alleleSeq_cnvScript and $ROOTSTOCK/allele/m116/m116_alleleSeq_cnvScript:
alleleWrap_2a_cnvnator_rd.sh
print_addRDcpp.sh
addCNV.cpp
Makefile


Run from cnvScript folder for given rootstock
NOTE: This conversion process is horribly slow...

\* well that's according to the manual, but addCNV (which creates the cnv format file) accepts the snv as well as the root file as input. I'm guessing it does some extra <i>stuff</i> which is taking the time.

```shell
$ROOTSTOCK/alleleSeq_cnvScript/alleleWrap_2a_cnvnator_rd.sh m27 $ROOTSTOCK/rootstock_genetics/ref/v1/fastas m27.root m27.snv

$ROOTSTOCK/alleleSeq_cnvScript/alleleWrap_2a_cnvnator_rd.sh m116 $ROOTSTOCK/rootstock_genetics/ref/v1/fastas m116.root m116.snv

```

### AlleleSeq pipeline

The following files were modified:
PIPELINE.mk (this is specific to each run)
MergeBowtie.py (set mappers to False as not looking for indels ~ this can probably be configured in PIPELINE.mk)
SnpCounts.py (set pat_1000G and mat_1000G to None for mappers and added line 48 chrm=chrm.replace("chr","") )

Updated allelseq pipeline to use Bowtie2 and accept paired-end read data. filter_input was modified as attached.

Trim RNA_Seq data

```shell
for d in $ROOTSTOCK/rna-seq/RNAseq*
do
 for f in ${d}/fastq/*1.fq.gz
 do
  r=$(echo $f|sed 's/_L1_1/_L1_2/')
  $ROOTSTOCK/scripts/fastq-mcf.sh $f $r $ROOTSTOCK/rootstock_genetics/ref/v1/Md ${d}/
 done
done
 ```
Created soft links to RNA-seq files in Alleleseq pipeline folder (links.sh).
Ran Alleleseq pipeline with following:

```shell
#rootstocks
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m27 03-RNA_L1_1.fq.trim 03-RNA_L1_2.fq.trim m27_paternal_index m27_maternal_index m27.snv m27.snv.cnv m27.map m27_S3_16B_hits.txt
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m27 05-RNA_L1_1.fq.trim 05-RNA_L1_2.fq.trim m27_paternal_index m27_maternal_index m27.snv m27.snv.cnv m27.map m27_S5_17B_hits.txt
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m27 07-RNA_L1_1.fq.trim 07-RNA_L1_2.fq.trim m27_paternal_index m27_maternal_index m27.snv m27.snv.cnv m27.map m27_S7_19B_hits.txt

$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m116 09-RNA_L1_1.fq.trim 09-RNA_L1_2.fq.trim m116_paternal_index m116_maternal_index m116.snv m116.snv.cnv m116.map m116_S9_6B_hits.txt
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m116 11-RNA_L1_1.fq.trim 11-RNA_L1_2.fq.trim m116_paternal_index m116_maternal_index m116.snv m116.snv.cnv m116.map m116_S11_8B_hits.txt
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m116 13-RNA_L1_1.fq.trim 13-RNA_L1_2.fq.trim m116_paternal_index m116_maternal_index m116.snv m116.snv.cnv m116.map m116_S13_10B_hits.txt

#gala
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m27 02-RNA_L1_1.fq.trim 02-RNA_L1_2.fq.trim m27_paternal_index m27_maternal_index m27.snv m27.snv.cnv m27.map m27_S2_16T_hits.txt
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m27 04-RNA_L1_1.fq.trim 04-RNA_L1_2.fq.trim m27_paternal_index m27_maternal_index m27.snv m27.snv.cnv m27.map m27_S4_17T_hits.txt
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m27 06-RNA_L1_1.fq.trim 06-RNA_L1_2.fq.trim m27_paternal_index m27_maternal_index m27.snv m27.snv.cnv m27.map m27_S6_19T_hits.txt

$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m116 08-RNA_L1_1.fq.trim 08-RNA_L1_2.fq.trim m116_paternal_index m116_maternal_index m116.snv m116.snv.cnv m116.map m116_S8_6T_hits.txt
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m116 10-RNA_L1_1.fq.trim 10-RNA_L1_2.fq.trim m116_paternal_index m116_maternal_index m116.snv m116.snv.cnv m116.map m116_S10_8T_hits.txt
$ROOTSTOCK/scripts/pipeline.sh $ROOTSTOCK/allele/m116 12-RNA_L1.1.fq.trim 12-RNA_L1.2.fq.trim m116_paternal_index m116_maternal_index m116.snv m116.snv.cnv m116.map m116_S12_10T_hits.txt
```



##QTL filtering


##Chromosome painting
http://www.paintmychromosomes.com/

Download and installed FineSTRUCTURE v2
```shell
tar -zxf fs-2.0.7.tar.gz
cd fs-2.0.7
./configure --prefix=$HOME/usr/local
make
make install
```
Beagle v4 output (vcf) is not compatible with FineSTRUCTURE and no conversion scripts are provided. Fortunately vcf2snp (see below) produces output which is easily converted to FineSTRUCTURE format.
FineSTRUCTURE format is as below

1. number samples (x2 for diploids)
2. number of snps
3. P followed by snp positions (space separted)
4. sample snp nucleotides maternal
5. sample snp nucleotides paternal
6. ++

e.g.  
4  
5  
P 1 5 10 20 30  
ATCGG  
ACCGT  
ATCGC  
CTCGG

The below was run to convert all snps in vcf2snp output files to FineSTRUCTURE format 
```shell
echo 10 > full.phase
wc -l m13.snv >>full.phase
echo -n "P " >>full.phase
cat mm106.snv|awk -F"\t" '{printf $2" "}'>>full.phase
echo >>full.phase
cat mm106.snv|awk -F"\t" '{printf substr ($6,1,1)}' >>full.phase
echo >>full.phase
cat mm106.snv|awk -F"\t" '{printf substr ($6,2,1)}' >>full.phase
echo >>full.phase
cat m13.snv|awk -F"\t" '{printf substr ($6,1,1)}' >>full.phase
echo >>full.phase
cat m13.snv|awk -F"\t" '{printf substr ($6,2,1)}' >>full.phase
echo >>full.phase
cat  m27.snv|awk -F"\t" '{printf substr ($6,1,1)}' >>full.phase
echo >>full.phase
cat  m27.snv|awk -F"\t" '{printf substr ($6,2,1)}' >>full.phase
echo >>full.phase
cat  m116.snv|awk -F"\t" '{printf substr ($6,1,1)}' >>full.phase
echo >>full.phase
cat  m116.snv|awk -F"\t" '{printf substr ($6,2,1)}' >>full.phase
echo >>full.phase
cat  m9.snv|awk -F"\t" '{printf substr ($6,1,1)}' >>full.phase
echo >>full.phase
cat m9.snv|awk -F"\t" '{printf substr ($6,2,1)}' >>full.phase
echo >>full.phase
```
To return snps from a subsection of one chromosome ( in this instance chr 5 between 1.2M and 36.5M bases) the following was run (using calc.pl):

```shell
echo 10 > chr5.int.phase
grep "^5" <m13.snv|perl calc.pl - 1200000 36500000|wc -l >>chr5.int.phase
echo -n "P " >>chr5.int.phase
grep "^5" mm106.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf $2" "}'>>chr5.int.phase
echo >>chr5.int.phase
grep "^5" mm106.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,1,1)}' >>chr5.int.phase
echo >>chr5.int.phase
grep "^5" mm106.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,2,1)}' >>chr5.int.phase
echo >>chr5.int.phase
grep "^5" m13.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,1,1)}' >>chr5.int.phase
echo >>chr5.int.phase
grep "^5" m13.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,2,1)}' >>chr5.int.phase
echo >>chr5.int.phase
grep "^5" m27.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,1,1)}' >>chr5.int.phase
echo >>chr5.int.phase
grep "^5" m27.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,2,1)}' >>chr5.int.phase
echo >>chr5.int.phase
grep "^5" m116.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,1,1)}' >>chr5.int.phase
echo >>chr5.int.phase
grep "^5" m116.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,2,1)}' >>chr5.int.phase
echo >>chr5.int.phase
grep "^5" m9.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,1,1)}' >>chr5.int.phase
echo >>chr5.int.phase
grep "^5" m9.snv|perl calc.pl - 1200000 36500000|awk -F"\t" '{printf substr ($6,2,1)}' >>chr5.int.phase
echo >>chr5.int.phase
```
To run FineSTRUCTURE

fs ch5.int.cp -idfile data.ids -phasefiles chr5.int.phase -go

(data.ids is a line seprated list of samples)

##Transcriptome ASE
http://www.cureffi.org/2013/08/26/allele-specific-rna-seq-pipeline-using-gsnap-and-gatk/
http://alleleseq.gersteinlab.org/downloads.html

The Allelseq pipeline, which uses vcf2diploid, then vcf2snp along with CNVnator 
looks like a good pipline to follow. 
The url above (gerstein) has most of the details and the readme files in vcf2snp is pretty comprehensive
##OLD COMMANDS


Made a shell script to submit automatically to grid engine
The reference genome was downloaded from GDR on 17/1/15
Reference location : /home/groups/harrisonlab/ref_genomes/rosaceae/malus/md_v2/Malus_x_domestica.v2.0-pht_assembly.fa

A Bowtie 2 hash was created

```shell
cd /home/groups/harrisonlab/project_files/rootstock_genetics/
mkdir ref
cd ref
bowtie2-build /home/groups/harrisonlab/ref_genomes/rosaceae/malus/md_v2/Malus_x_domestica.v2.0-pht_assembly.fa Md
```

```shell
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m27/conc/m27_r1.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/m27/conc/m27_r2.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m27/analysis/ m27_ge.sam 250 500
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m116/conc/m116_r1.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/m116/conc/m116_r2.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m116/analysis/ m116_ge.sam 250 500
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m9/conc/m9_r1.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/m9/conc/m9_r2.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m9/analysis/ m9_ge.sam 250 500
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m13/conc/m13_r1.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/m13/conc/m13_r2.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m13/analysis/ m13_ge.sam 250 500 
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/conc/mm106_r1.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/conc/mm106_r2.fq.trim /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/analysis/ mm106_ge.sam 250 500
```
Convert SAM to BAM and sort BAM for SNP calling
Do this by SGE script
```shell
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m27/analysis/m27_ge.sam /home/groups/harrisonlab/project_files/rootstock_genetics/m27/analysis/ m27.bam m27.sorted 
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m9/analysis/m9_ge.sam /home/groups/harrisonlab/project_files/rootstock_genetics/m9/analysis/ m9.bam m9.sorted
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m116/analysis/m116_ge.sam /home/groups/harrisonlab/project_files/rootstock_genetics/m116/analysis/ m116.bam m116.sorted
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m13/analysis/m13_ge.sam /home/groups/harrisonlab/project_files/rootstock_genetics/m13/analysis/ m13.bam m13.sorted 
./samtools.sh /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/analysis/mm106_ge.sam /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/analysis/ mm106.bam mm106.sorted

```
Index the reference genome for SAMTOOLS
```shell
samtools faidx Malus_x_domestica.v2.0-pht_assembly.fa 

```
Pileup into a single vcf

```shell
#samtools mpileup -uf ./ref/Malus_x_domestica.v2.0-pht_assembly.fa.fai  ./m9/analysis/m9-sorted.bam  ./m27/analysis/m27-sorted.bam  ./m116/analysis/m116-sorted.bam  ./mm106/analysis/mm106-sorted.bam  ./m13/analysis/m13-sorted.bam | bcftools view -bvcg - > ./vcf/var.raw.bcf

samtools mpileup -uf ./ref/Malus_x_domestica.v2.0-pht_assembly.fa  ./m9/analysis/m9.sorted.bam  ./m27/analysis/m27.sorted.bam  ./m116/analysis/m116.sorted.bam  ./mm106/analysis/mm106.sorted.bam  ./m13/analysis/m13.sorted.bam >pileup.bam
bcftools view -bvcg pileup.bam > ./vcf/var.raw.bcf
bcftools view ./vcf/var.raw.bcf | vcfutils.pl varFilter -D100 > ./beagle/var.flt.vcf

```


Trim ASE 
```shell
./fastq-mcf.sh $ROOTSTOCK/rna-seq/RNAseq_3/02-RNA_L1_1.fq $ROOTSTOCK/rna-seq/RNAseq_3/02-RNA_L1_2.fq $ROOTSTOCK/rootstock_genetics/ref/v1/Md $ROOTSTOCK/rna-seq/RNAseq_2/ 
./fastq-mcf.sh $ROOTSTOCK/rna-seq/RNAseq_3/03-RNA_L1_1.fq $ROOTSTOCK/rna-seq/RNAseq_3/03-RNA_L1_2.fq $ROOTSTOCK/rootstock_genetics/ref/v1/Md $ROOTSTOCK/rna-seq/RNAseq_2/ 
./fastq-mcf.sh $ROOTSTOCK/rna-seq/RNAseq_3/04-RNA_L1_1.fq $ROOTSTOCK/rna-seq/RNAseq_3/04-RNA_L1_2.fq $ROOTSTOCK/rootstock_genetics/ref/v1/Md $ROOTSTOCK/rna-seq/RNAseq_2/ 
./fastq-mcf.sh $ROOTSTOCK/rna-seq/RNAseq_3/05-RNA_L1_1.fq $ROOTSTOCK/rna-seq/RNAseq_3/05-RNA_L1_2.fq $ROOTSTOCK/rootstock_genetics/ref/v1/Md $ROOTSTOCK/rna-seq/RNAseq_2/
./fastq-mcf.sh $ROOTSTOCK/rna-seq/RNAseq_2/06-RNA_L1_1.fq $ROOTSTOCK/rna-seq/RNAseq_2/06-RNA_L1_2.fq $ROOTSTOCK/rootstock_genetics/ref/v1/Md $ROOTSTOCK/rna-seq/RNAseq_2/ 
./fastq-mcf.sh $ROOTSTOCK/rna-seq/RNAseq_2/07-RNA_L1_1.fq $ROOTSTOCK/rna-seq/RNAseq_2/07-RNA_L1_2.fq $ROOTSTOCK/rootstock_genetics/ref/v1/Md /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/ 
./fastq-mcf.sh /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/08-RNA_L1_1.fq /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/08-RNA_L1_2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/ 
./fastq-mcf.sh /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/09-RNA_L1_1.fq /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/09-RNA_L1_2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/ 
./fastq-mcf.sh /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_1/10-RNA_L1_1.fq /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_1/10-RNA_L1_2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/ 
./fastq-mcf.sh /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_1/11-RNA_L1.1.fq /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_1/11-RNA_L1.2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/ 
./fastq-mcf.sh /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_1/12-RNA_L1.1.fq /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_1/12-RNA_L1.2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/ 
./fastq-mcf.sh /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_1/12-RNA_L1.1.fq /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_1/13-RNA_L1.2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/v1/Md /home/deakig/projects/apple_rootstock/rna-seq/RNAseq_2/ 
```
