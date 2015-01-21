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

nohup bowtie2 -p 4 -x ./ref/Md -1 ./m27/conc/m27_r1.fq -2 ./m27/conc/m27_r2.fq -S ./m27/analysis/m27.sam &>m27.out&
nohup bowtie2 -p 4 -x ./ref/Md -1 ./m116/conc/m116_r1.fq -2 ./m116/conc/m116_r2.fq -S ./m116/analysis/m116.sam &>m116.out&
nohup bowtie2 -p 4 -x ./ref/Md -1 ./m9/conc/m9_r1.fq -2 ./m9/conc/m9_r2.fq -S ./m9/analysis/m9.sam &>m9.out&
nohup bowtie2 -p 4 -x ./ref/Md -1 ./m13/conc/m13_r1.fq -2 ./m13/conc/m13_r2.fq -S ./m13/analysis/m13.sam &>m13.out&
nohup bowtie2 -p 4 -x ./ref/Md -1 ./mm106/conc/mm106_r1.fq -2 ./mm106/conc/mm106_r2.fq -S ./mm106/analysis/mm106.sam &>mm106.out&

```

Made a shell script to submit automatically to grid engine

```shell
./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m27/conc/m27_r1.fq /home/groups/harrisonlab/project_files/rootstock_genetics/m27/conc/m27_r2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m27/analysis/m27_ge.sam

./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m116/conc/m116_r1.fq /home/groups/harrisonlab/project_files/rootstock_genetics/m116/conc/m116_r2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m116/analysis/m116_ge.sam

./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m9/conc/m9_r1.fq /home/groups/harrisonlab/project_files/rootstock_genetics/m9/conc/m9_r2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m9/analysis/m9_ge.sam

./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/m13/conc/m13_r1.fq /home/groups/harrisonlab/project_files/rootstock_genetics/m13/conc/m13_r2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/m13/analysis/m13_ge.sam

./bowtie.sh /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/conc/mm106_r1.fq /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/conc/mm106_r2.fq /home/groups/harrisonlab/project_files/rootstock_genetics/ref/Md /home/groups/harrisonlab/project_files/rootstock_genetics/mm106/analysis/mm106_ge.sam
```



Convert SAM to BAM for sorting
```shell
nohup samtools view -S -b ./m9/analysis/m9_ge.sam >./m9/analysis/m9.bam &
nohup samtools view -S -b ./m27/analysis/m27_ge.sam >./m27/analysis/m27.bam  &
nohup samtools view -S -b ./m116/analysis/m116_ge.sam >./m116/analysis/m116.bam  &
nohup samtools view -S -b ./mm106/analysis/mm106_ge.sam >./mm106/analysis/mm106.bam  &
nohup samtools view -S -b ./m13/analysis/m13_ge.sam >./m13/analysis/m13.bam  &

```
 Sort BAM for SNP calling
```shell
nohup samtools sort ./m9/analysis/m9.bam ./m9/analysis/m9-sorted.bam  &
nohup samtools sort ./m27/analysis/m27.bam ./m27/analysis/m27-sorted.bam &
nohup samtools sort ./m116/analysis/m116.bam ./m116/analysis/m116-sorted.bam &
nohup samtools sort ./mm106/analysis/mm106.bam ./mm106/analysis/mm106-sorted.bam &
nohup samtools sort ./m13/analysis/m13.bam ./m13/analysis/m13-sorted.bam &
```
Index the reference genome for SAMTOOLS
```shell
samtools faidx Malus_x_domestica.v2.0-pht_assembly.fa 
```

Pileup into a single vcf

```shell
samtools mpileup -uf Malus_x_domestica.v2.0-pht_assembly.fa.fai  m9-sorted.bam  m27-sorted.bam  m116-sorted.bam  mm106-sorted.bam  m13-sorted.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > ./beagle/var.flt.vcf
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

java –Xmx 12000m –jar beagle.jar gt=./beagle/var.flt.vcf ped=./beagle/pedigree.ped out=beagle chrom=1 nthreads=4


samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf


http://www.paintmychromosomes.com/

##Transcriptome ASE
http://www.cureffi.org/2013/08/26/allele-specific-rna-seq-pipeline-using-gsnap-and-gatk/
http://alleleseq.gersteinlab.org/downloads.html

The Allelseq pipeline, which uses vcf2diploid, then vcf2snp along with CNVnator 
looks like a good pipline to follow. 
The url above (gerstein) has most of the details and the readme files in vcf2snp is pretty comprehensive


##QTL filtering



