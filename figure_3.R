options(stringsAsFactors = FALSE)
library(plyr)

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
	#write.table(unique_phased_genes,xx,sep="\t",quote=F,row.names=F)
	mylist <- list(t1,unique_phased_genes)
	names(mylist) <- c("allele_seq","sig_phased")
	return(mylist)
}

testfunc <- function(x,Y) {
# determines the proportion of allelic expressed snps within a gene which are on the maternal chromosome
	Y <- Y[Y$pos>=x[3]&Y$pos<=x[4]&Y$evidence>=2,3:8]
	colnames(Y)[3:6] <- c("A","C","G","T")
	if(!is.null(Y)) {
		test <- (sum(Y[,1]==colnames(Y[3:6])[max.col(Y[3:6],ties.method="random")])/length(Y[,1]))
		return(round(test,2))
	} else {
		return(NA)
	}
}

exon_correction <- function(x,Y) {
# get_data sums the no. of snps per gene - this is a quick fix to get the per exon values
	Y <- Y[Y$pos>=x[3]&Y$pos<=x[4]&Y$evidence>=2,3:8]
	return(length(Y[,1]))
}

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

# ~/projects/apple_rootstock/allele/Figure3/hits
# m27_chr5.txt
# m27_chr11.txt
# m27_chr13.txt
# m116_chr5.txt - dw1 is not present in m116
# m116_chr11.txt
# m116_chr13.txt
# abovoe files filtered and corrected with correct_phase.r to give files with the same name minus the suffix

# get data for genes
m27_pg_q5 <-  get_data("m27_chr5","qtl5_genes",0.05)
m27_pg_q11 <- get_data("m27_chr11","qtl11_genes",0.05)
m27_pg_q13 <- get_data("m27_chr13","qtl13_genes",0.05)
m116_pg_q5 <- get_data("m116_chr5","qtl5_genes",0.05)
m116_pg_q11 <- get_data("m116_chr11","qtl11_genes",0.05)
m116_pg_q13 <- get_data("m116_chr13","qtl13_genes",0.05)

# get the proporion of allelic expressed snps (per gene) which are on the maternal chromosome
m27_pg_q5$sig_phased$mat_prop <- apply(m27_pg_q5$sig_phased,1,function(x) testfunc(x,m27_pg_q5$allele_seq))
m27_pg_q11$sig_phased$mat_prop <- apply(m27_pg_q11$sig_phased,1,function(x) testfunc(x,m27_pg_q11$allele_seq))
m27_pg_q13$sig_phased$mat_prop <- apply(m27_pg_q13$sig_phased,1,function(x) testfunc(x,m27_pg_q13$allele_seq))
m116_pg_q5$sig_phased$mat_prop <- apply(m116_pg_q5$sig_phased,1,function(x) testfunc(x,m116_pg_q5$allele_seq))
m116_pg_q11$sig_phased$mat_prop <- apply(m116_pg_q11$sig_phased,1,function(x) testfunc(x,m116_pg_q11$allele_seq))
m116_pg_q13$sig_phased$mat_prop <- apply(m116_pg_q13$sig_phased,1,function(x) testfunc(x,m116_pg_q13$allele_seq))

# get data for exons
m27_pe_q5 <- get_data("m27_chr5","qtl5_exons",0.05)
m27_pe_q11 <- get_data("m27_chr11","qtl11_exons",0.05)
m27_pe_q13 <- get_data("m27_chr13","qtl13_exons",0.05)
m116_pe_q5 <- get_data("m116_chr5","qtl5_exons",0.05)
m116_pe_q11 <- get_data("m116_chr11","qtl11_exons",0.05)
m116_pe_q13 <- get_data("m116_chr13","qtl13_exons",0.05)
# oops - get_data sums the number of snps across genes - exon_correction will correct this
m27_pe_q5$sig_phased$snp_count <- apply(m27_pe_q5$sig_phased,1,function(x) exon_correction(x,m27_pe_q5$allele_seq))
m27_pe_q11$sig_phased$snp_count <- apply(m27_pe_q11$sig_phased,1,function(x) exon_correction(x,m27_pe_q11$allele_seq))
m27_pe_q13$sig_phased$snp_count <- apply(m27_pe_q13$sig_phased,1,function(x) exon_correction(x,m27_pe_q13$allele_seq))
m116_pe_q5$sig_phased$snp_count <- apply(m116_pe_q5$sig_phased,1,function(x) exon_correction(x,m116_pe_q5$allele_seq))
m116_pe_q11$sig_phased$snp_count <- apply(m116_pe_q11$sig_phased,1,function(x) exon_correction(x,m116_pe_q11$allele_seq))
m116_pe_q13$sig_phased$snp_count <- apply(m116_pe_q13$sig_phased,1,function(x) exon_correction(x,m116_pe_q13$allele_seq))
# Remove exons without allelic snps
m27_pe_q5$sig_phased <- m27_pe_q5$sig_phased[m27_pe_q5$sig_phased$snp_count>0,]
m27_pe_q11$sig_phased <- m27_pe_q11$sig_phased[m27_pe_q11$sig_phased$snp_count>0,]
m27_pe_q13$sig_phased <- m27_pe_q13$sig_phased[m27_pe_q13$sig_phased$snp_count>0,]
m116_pe_q5$sig_phased <- m116_pe_q5$sig_phased[m116_pe_q5$sig_phased$snp_count>0,]
m116_pe_q11$sig_phased <- m116_pe_q11$sig_phased[m116_pe_q11$sig_phased$snp_count>0,]
m116_pe_q13$sig_phased <- m116_pe_q13$sig_phased[m116_pe_q13$sig_phased$snp_count>0,]
# get the proporion of allelic expressed snps (per exon) which are on the maternal chromosome
m27_pe_q5$sig_phased$mat_prop <- apply(m27_pe_q5$sig_phased[,1:4],1,function(x) testfunc(x,m27_pe_q5$allele_seq))
m27_pe_q11$sig_phased$mat_prop <- apply(m27_pe_q11$sig_phased[,1:4],1,function(x) testfunc(x,m27_pe_q11$allele_seq))
m27_pe_q13$sig_phased$mat_prop <- apply(m27_pe_q13$sig_phased[,1:4],1,function(x) testfunc(x,m27_pe_q13$allele_seq))
m116_pe_q5$sig_phased$mat_prop <- apply(m116_pe_q5$sig_phased[,1:4],1,function(x) testfunc(x,m116_pe_q5$allele_seq))
m116_pe_q11$sig_phased$mat_prop <- apply(m116_pe_q11$sig_phased[,1:4],1,function(x) testfunc(x,m116_pe_q11$allele_seq))
m116_pe_q13$sig_phased$mat_prop <- apply(m116_pe_q13$sig_phased[,1:4],1,function(x) testfunc(x,m116_pe_q13$allele_seq))

# combine gene and exon results
combine <- function(x,Y) {
	Y <- Y[Y$Gene_ID==x[1],]
	Y <- Y[order(Y$Start),]
	return(paste(Y[,5],":",Y[,2],":",Y[,7],sep="",collapse=";"))
}

m27_pg_q5$sig_phased$exons <- apply(m27_pg_q5$sig_phased,1,function(x) combine(x,m27_pe_q5$sig_phased))
m27_pg_q11$sig_phased$exons <- apply(m27_pg_q11$sig_phased,1,function(x) combine(x,m27_pe_q11$sig_phased))
m27_pg_q13$sig_phased$exons <- apply(m27_pg_q13$sig_phased,1,function(x) combine(x,m27_pe_q13$sig_phased))
m116_pg_q5$sig_phased$exons <- apply(m116_pg_q5$sig_phased,1,function(x) combine(x,m116_pe_q5$sig_phased))
m116_pg_q11$sig_phased$exons <- apply(m116_pg_q11$sig_phased,1,function(x) combine(x,m116_pe_q11$sig_phased))
m116_pg_q13$sig_phased$exons <- apply(m116_pg_q13$sig_phased,1,function(x) combine(x,m116_pe_q13$sig_phased))

# write output
write.table(m27_pg_q5$sig_phased,"m27_q5_allelic_genes",sep="\t",row.names=F,quote=F)
write.table(m27_pg_q11$sig_phased,"m27_q11_allelic_genes",sep="\t",row.names=F,quote=F)
write.table(m27_pg_q13$sig_phased,"m27_q13_allelic_genes",sep="\t",row.names=F,quote=F)
write.table(m116_pg_q5$sig_phased,"m116_q5_allelic_genes",sep="\t",row.names=F,quote=F)
write.table(m116_pg_q11$sig_phased,"m116_q11_allelic_genes",sep="\t",row.names=F,quote=F)
write.table(m116_pg_q13$sig_phased,"m116_q13_allelic_genes",sep="\t",row.names=F,quote=F)

##  stuff
q5_uncom <-  m27_pg_q5$sig_phased[!m27_pg_q5$sig_phased[,1]%in%m116_pg_q5$sig_phased[,1],] 
q11_com <- merge(m27_pg_q11$sig_phased,m116_pg_q11$sig_phased,by="Gene_ID")
q13_com <- merge(m27_pg_q13$sig_phased,m116_pg_q13$sig_phased,by="Gene_ID")


