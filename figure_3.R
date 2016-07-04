options(stringsAsFactors = FALSE)
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

ev <- function(p,a) {
	if(p<=a) {
		return(1)
	} else {
		return(0)
	}
}

genes <- read.table("qtl13_genes",header=T,sep="\t")
exons <- read.table("qtl_exons",header=T,sep="\t")

#m27_chr5.txt
#m27_chr11.txt
#m27_chr13.txt
#m116_chr5.txt - dw1 is not present in m116
#m116_chr11.txt
#m116_chr13.txt

alpha = 0.05
x="m27_chr13"
t1 <- read.table(paste(x,".txt",sep=""),header=T,sep="\t")

t1$evidence <- apply(t1[,12:14],1,function(x) return(ev(x[1],alpha)+ev(x[2],alpha)+ev(x[3],alpha)))
phased <- t1[t1[,2]=="PHASED"&t1$evidence>=2,] 

#list of genes and exons containing phased snps

##produces list of data.frames of 0 or greater rows
phased_genes <- apply(phased,1,function(x) within(genes,x[1]))
##combine into a single data.frame (dumps 0 length rows)
phased_genes <- ldply(phased_genes,data.frame)
##collapse to  unique plus number of snps per gene
unique_phased_genes <- ddply(phased_genes,.(Gene_ID),nrow)
colnames(unique_phased_genes)[2] <- "snp_count"
unique_phased_genes <- merge(unique_phased_genes,genes,all.x=T)
xx <- paste(x,"_phased_genes.txt",sep="")
write.table(unique_phased_genes,xx,sep="\t",quote=F,row.names=F)

##get unique phased exons
phased_exons <- apply(phased,1,function(x) within(exons,x[1]))
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


