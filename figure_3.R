options(stringsAsFactors = FALSE)
library(plyr)

between <- function(X,pos,dist) {
	X[(X[,1]>=(pos-dist)&X[,2]<=(pos+dist)),]
}

within <- function(X,pos,bin) {
	X[(X[,1]<=pos&X[,2]>=pos),]
}

upstream <- function(X,pos,dist) {
	X[
		(X$direction=="forward"&pos<X[,1]&pos>=(X[,1]-dist)) |
		(X$direction=="reverse"&pos>X[,2]&pos<=(X[,2]+dist))
	,]
}

downstream <- function(X,pos,dist) {
	X[
		(X$direction=="forward"&pos>X[,2]&pos<=(X[,2]+dist)) |
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

get_data <- function(sname,loc,alpha=0.05,location=within,dist=500,evidence=2) {
	options(stringsAsFactors = FALSE)
	library(plyr)
	mytable <- read.table(sname,header=T,sep="\t")
	#genes could also be a list of exons
	genes <- read.table(loc,header=T,sep="\t")
	# Set probability of Weird results to 1 (they're too weird to use)
	t1 <- mytable
	t1[grep("Weird",t1[,9]),12] <- 1
	t1[grep("Weird",t1[,10]),13] <- 1
	t1[grep("Weird",t1[,11]),14] <- 1
	# sum samples with allele-seq probability less than or equal to alpha 
	#alpha = 0.05
	t1$evidence <- apply(t1[,12:14],1,function(x) return(ev(x[1],alpha)+ev(x[2],alpha)+ev(x[3],alpha)))
	# Select phased snps with at least 2 samples with significant allelic expression 
	# (selecting on column 2 (phase) is not strictly necessary as only phased snps will have a probability less than 1)
	phased <- t1[t1[,2]=="PHASED"&t1$evidence>=evidence,] 
	#list of genes and exons containing phased snps
	##produces list of data.frames of 0 or greater rows
	phased_genes <- apply(phased,1,function(x) location(genes,x[1],dist))
	##combine into a single data.frame (dumps 0 length rows)
	phased_genes <- ldply(phased_genes,data.frame)
	##collapse to  unique plus number of snps per gene
	unique_phased_genes <- ddply(phased_genes,.(Gene_ID),nrow)
	colnames(unique_phased_genes)[2] <- "snp_count"
	unique_phased_genes <- merge(unique_phased_genes,genes,all.x=T)
	xx <- paste(sname,"_phased_genes.txt",sep="")	
	write.table(unique_phased_genes,xx,sep="\t",quote=F,row.names=F)
	mylist <- list(t1,unique_phased_genes)
	names(mylist) <- c("allele_seq","sig_phased")
	return(mylist)
}

exons <- read.table("qtl_exons",header=T,sep="\t")

#m27_chr5.txt
#m27_chr11.txt
#m27_chr13.txt
#m116_chr5.txt - dw1 is not present in m116
#m116_chr11.txt
#m116_chr13.txt

m27_pg_q5 <-  get_data("m27_chr5","qtl5_genes",0.05)
m27_pg_q11 <- get_data("m27_chr11","qtl11_genes",0.05)
m27_pg_q13 <- get_data("m27_chr13","qtl13_genes",0.05)
m116_pg_q5 <- get_data("m116_chr5","qtl5_genes",0.05)
m116_pg_q11 <- get_data("m116_chr11","qtl11_genes",0.05)
m116_pg_q13 <- get_data("m116_chr13","qtl13_genes",0.05)

test <- apply(m27_pg_q5$sig_phased,1,function(x) x)
gitfunc <- function(x,Y) {
	Y <- Y[Y$pos>=x[3]&Y$pos<=x[4]&Y&evidence>=2,3:8]
	colnames(Y)[3:6] <- C("A","C","G","T")

	Y[3],Y[4] vs Y[5][6][7][8]
}


test <- m27_pg_q5$allele_seq[m27_pg_q5$allele_seq$pos>=m27_pg_q5$sig_phased$start[1]&m27_pg_q5$allele_seq$pos<=m27_pg_q5$sig_phased$end[1],]

q5_uncom <-  m27_pg_q5$sig_phased[!m27_pg_q5$sig_phased[,1]%in%m116_pg_q5$sig_phased[,1],] 
q11_com <- merge(m27_pg_q11$sig_phased,m116_pg_q11$sig_phased,by="Gene_ID")
q13_com <- merge(m27_pg_q13$sig_phased,m116_pg_q13$sig_phased,by="Gene_ID")

##get unique phased exons
#phased_exons <- apply(phased,1,function(x) within(exons,x[1]))
#phased_exons <- ldply(phased_exons ,data.frame)

##collapse to  unique plus number of snps per gene (i.e. combine counts from exons)
#unique_phased_exons <- ddply(phased_exons,.(Gene_ID),nrow)
#colnames(unique_phased_exons)[2] <- "snp_count"

# genes and exons without phased snps??

###gene_snps <- apply(snps_phased[1:100,],1,function(x) within(genes,x[2])) --- way too slow

#gene_snps <- apply(genes,1,function(x) as.data.frame(cbind(x[3],sum(as.logical(snps_phased$V2>=as.integer(x[1])&snps_phased$V2<=as.integer(x[2]))))))
#gene_snps <- ldply(gene_snps,data.frame)
#dim(gene_snps[as.integer(levels(gene_snps$V2))[gene_snps$V2]>0,])



