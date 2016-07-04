library(plyr)


between <- function(X,pos,dist) {
	X[
		(X$direction=="forward"&X[,1]>=(pos-dist)&X[,2]<=(pos+dist)) |
		(X$direction=="reverse"&X[,2]>=(pos-dist)&X[,1]<=(pos+dist))
	,]
}

within <- function(X,pos) {
	X[
		(X$direction=="forward"&X[,1]<=pos&X[,2]>=pos) |
		(X$direction=="reverse"&X[,2]<=pos&X[,1]>=pos)
	,]
}

upstream <- function(X,pos,dist) {
	X[
		(X$direction=="forward"&pos<X[,1]&pos>=(X[,1]-dist)) |
		(X$direction=="reverse"&pos>X[,2]&pos<=(X[,2]-dist))
	,]
}

downstream <- function(X,pos,dist) {
	X[
		(X$direction=="forward"&pos>X[,2]&pos<=(X[,2]-dist)) |
		(X$direction=="reverse"&pos<X[,1]&pos>=(X[,1]-dist))
	,]
}


genes <- read.table("qtl_genes",header=T,sep="\t")
exons <- read.table("qtl_exons",header=T,sep="\t")
phased <- read.table("phased",header=T,sep="\t") ## phased snps showing significant allelic expression differences
snps <- read.table("../snps/ver2/m27.snv",header=F,sep="\t")
snps <- snps[snps$V1==5,]


#list of genes and exons containing phased snps

##produces list of data.frames of 0 or greater rows
phased_genes <- apply(phased,1,function(x) within(genes,x[2]))
##combine into a single data.frame (dumps 0 length rows)
phased_genes <- ldply(phased_genes,data.frame)
##collapse to  unique plus number of snps per gene
unique_phased_genes <- ddply(phased_genes,.(Gene_ID),nrow)
colnames(unique_phased_genes)[2] <- "snp_count"

##get unique phased exons
phased_exons <- apply(phased,1,function(x) within(exons,x[2]))
phased_exons <- ldply(phased_exons ,data.frame)

##collapse to  unique plus number of snps per gene (i.e. combine counts from exons)
unique_phased_exons <- ddply(phased_exons,.(Gene_ID),nrow)
colnames(unique_phased_exons)[2] <- "snp_count"

# genes and exons without phased snps??

###gene_snps <- apply(snps_phased[1:100,],1,function(x) within(genes,x[2])) --- way too slow

gene_snps <- apply(genes,1,function(x) as.data.frame(cbind(x[3],sum(as.logical(snps_phased$V2>=as.integer(x[1])&snps_phased$V2<=as.integer(x[2]))))))
gene_snps <- ldply(gene_snps,data.frame)




#### Results QTL5 - m27
QTL_genes 		= 1786
dim(gene_snps[as.integer(levels(gene_snps$V2))[gene_snps$V2]>0,])
genes with phased snps 	= 1266
genes sig allele expressed	= 199 (130 within exons)


