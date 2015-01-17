# rootstock

Data was unzipped and concatonated into R1 and R2 files in the subdirectory conc within each rootstock folder
```shell
gunzip *
cat 863_LIB6292_LDI5172_GTGAAA_L00*_R1.fastq.gz >m27_r1.fa
``
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

To generate a whole consensus sequence
samtools mpileup -uf ref.fa aln.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq i 

##Phasing with Beagle
http://faculty.washington.edu/browning/beagle_utilities/utilities.html

samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf


##Transcriptome ASE


##QTL filtering



