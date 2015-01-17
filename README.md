# rootstock

To run this analysis I used:
fastqc
bowtie2
samtools
bcftools
beagle



Data was unzipped and concatonated into R1 and R2 files in the subdirectory conc within each rootstock folder
```shell
gunzip *
cat 863_LIB6292_LDI5172_GTGAAA_L00*_R1.fastq.gz >m27_r1.fa
```
fastqc was run on each pair of reads in order to assess the quality
```shell
nohup fastqc m116_r1.fq m116_r2.fq &
nohup fastqc m27_r1.fq m27_r2.fq &
nohup fastqc m9_r1.fq m9_r2.fq &
nohup fastqc m13_r1.fq m13_r2.fq &
nohup fastqc mm106_r1.fq mm106_r2.fq &
 ```
##Assembly to reference

The reference genome was downloaded from GDR on 17/1/15

Reference location : /home/groups/harrisonlab/ref_genomes/rosaceae/malus/md_v2/Malus_x_domestica.v2.0-pht_assembly.fa

A Bowtie 2 hash was created

```shell
cd /home/groups/harrisonlab/project_files/rootstock_genetics/
mkdir ref
cd ref
bowtie2-build /home/groups/harrisonlab/ref_genomes/rosaceae/malus/md_v2/Malus_x_domestica.v2.0-pht_assembly.fa Md

nohup bowtie2 -x ./ref/Md -1 ./m27/conc/m27_r1.fq -2 ./m27/conc/m27_r2.fq -S ./m27/analysis/m27.sam &>m27.out&
nohup bowtie2 -x ./ref/Md -1 ./m116/conc/m116_r1.fq -2 ./m116/conc/m116_r2.fq -S ./m116/analysis/m116.sam &>m116.out&
nohup bowtie2 -x ./ref/Md -1 ./m9/conc/m9_r1.fq -2 ./m9/conc/m9_r2.fq -S ./m9/analysis/m9.sam &>m9.out&
nohup bowtie2 -x ./ref/Md -1 ./m13/conc/m13_r1.fq -2 ./m13/conc/m13_r2.fq -S ./m13/analysis/m13.sam &>m13.out&
nohup bowtie2 -x ./ref/Md -1 ./mm106/conc/mm106_r1.fq -2 ./mm106/conc/mm106_r2.fq -S ./mm106/analysis/mm106.sam &>mm106.out&

```


Convert SAM to BAM for sorting
```shell
nohup samtools view -S -b ./m9/analysis/m9.sam >./m9/analysis/m9.bam &>m9.samtools.out&
nohup samtools view -S -b ./m27/analysis/m27.sam >./m27/analysis/m27.bam  &>m27.samtools.out&
nohup samtools view -S -b ./m116/analysis/m116.sam >./m116/analysis/m116.bam  &>m116.samtools.out&
nohup samtools view -S -b ./mm106/analysis/mm106.sam >./mm106/analysis/mm106.bam  &>mm106.samtools.out&
nohup samtools view -S -b ./m13/analysis/m13.sam >./m13/analysis/m13.bam  &>m13.samtools.out&

```
 Sort BAM for SNP calling
```shell
nohup samtools sort ./m9/analysis/m9.bam ./m9/analysis/m9-sorted.bam  &>m9.sorted.out&
nohup samtools sort ./m9/analysis/m27.bam ./m9/analysis/m27-sorted.bam &>m27.sorted.out&
nohup samtools sort ./m9/analysis/m116.bam ./m9/analysis/m116-sorted.bam &>m116.sorted.out&
nohup samtools sort ./m9/analysis/mm106.bam ./m9/analysis/mm106-sorted.bam &>mm106.sorted.out&
nohup samtools sort ./m9/analysis/m13.bam ./m9/analysis/m13-sorted.bam &>m13.sorted.out&
```
Index the reference genome for SAMTOOLS
```shell
samtools faidx Malus_x_domestica.v2.0-pht_assembly.fa 
```

Pileup into a single vcf

```shell
samtools mpileup -uf Malus_x_domestica.v2.0-pht_assembly.fa.fai  m9-sorted.bam  m27-sorted.bam  m116-sorted.bam  mm106-sorted.bam  m13-sorted.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf
```

Define the pedigree for beagle
1) pedigree ID, 2) individual ID, 3) father’s ID, and 4) mother’s ID

m116 1 2 3
mm106 2 0 0
m27 3 4 5
m13 4 0 0
m9 5 0 0


NB- To generate a whole consensus sequence
samtools mpileup -uf ref.fa aln.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq i 

##Phasing with Beagle
http://faculty.washington.edu/browning/beagle_utilities/utilities.html
Note SHAPE could also be used- one advantage here is it can be 'read aware'  https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#readaware


samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf


##Transcriptome ASE


##QTL filtering



