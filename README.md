# rootstock

##Assembly to reference

Reference: /home/groups/harrisonlab/ref_genomes/rosaceae/malus/md_v2/Malus_x_domestica.v2.0-pht_assembly.fa
```shell
cd /home/groups/harrisonlab/project_files/rootstock_genetics/
mkdir ref
cd ref
bowtie2-build /home/groups/harrisonlab/ref_genomes/rosaceae/malus/md_v2/Malus_x_domestica.v2.0-pht_assembly.fa Md

bowtie2 -x ./ref/Md -U ./M27/reads_1.fq -S m27.sam


samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf  
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf  
```
To generate a whole consensus sequence
samtools mpileup -uf ref.fa aln.bam | bcftools view -cg - | vcfutils.pl vcf2fq > cns.fq i 

##Phasing with Beagle
http://faculty.washington.edu/browning/beagle_utilities/utilities.html

##Transcriptome ASE


##QTL filtering



