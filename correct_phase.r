#The phasing in the allele-seq pipeline may not be ideal
#This script will merge the allele-seq output with the results from phasing.r


quickfunc <- function(v) {
  line <-  c(v[1],strsplit(v[2],","))
  return(line((as.numeric(v[3])+1)))
}

#read in output from stuff(fig_3).r
m27_chr5 <- read.table("m27_chr5.txt",header=T,sep="\t") 
#read in phased_m27 snps from phasing.r
phased_m27 <- read.table("phased_m27",sep="\t",header=T)
#merge data
phased_m27_chr5 <- merge(m27_chr5[m27_chr5$phase=="PHASED",],phased_m27[(phased_m27$chr=="chr5"),],by.x="pos",by.y="pos")
#set mat nucleotide to correct phase ( head(phased_m27_chr5) shows that the 5th snp is incorrect)
phased_m27_chr5$mat <- apply(phased_m27_chr5[,16:18],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
#set pat nucleotide to correct phase
phased_m27_chr5$pat <- apply(phased_m27_chr5[,c(16,17,19)],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
#remove unnecessary rows
phased_m27_chr5 <- phased_m27_chr5[,1:14]
#write output
write.table(phased_m27_chr5,"m27_chr5",sep="\t",row.names=F,quote=F)

# and the rest
Y <- phased_m27 
#m27chr11
X <- read.table("m27_chr11.txt",header=T,sep="\t") 
P <- merge(X[X$phase=="PHASED",],Y[(Y$chr=="chr11"),],by.x="pos",by.y="pos")
P$mat <- apply(P[,16:18],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P$pat <- apply(P[,c(16,17,19)],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P <- P[,1:14]
write.table(P,"m27_chr11",sep="\t",row.names=F,quote=F)
#m27chr13
X <- read.table("m27_chr13.txt",header=T,sep="\t") 
P <- merge(X[X$phase=="PHASED",],Y[(Y$chr=="chr11"),],by.x="pos",by.y="pos")
P$mat <- apply(P[,16:18],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P$pat <- apply(P[,c(16,17,19)],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P <- P[,1:14]
write.table(P,"m27_chr13",sep="\t",row.names=F,quote=F)

Y <- phased_m116 
#m116chr5
X <- read.table("m116_chr5.txt",header=T,sep="\t") 
P <- merge(X[X$phase=="PHASED",],Y[(Y$chr=="chr11"),],by.x="pos",by.y="pos")
P$mat <- apply(P[,16:18],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P$pat <- apply(P[,c(16,17,19)],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P <- P[,1:14]
write.table(P,"m116_chr5",sep="\t",row.names=F,quote=F)
#m116chr11
X <- read.table("m116_chr11.txt",header=T,sep="\t") 
P <- merge(X[X$phase=="PHASED",],Y[(Y$chr=="chr11"),],by.x="pos",by.y="pos")
P$mat <- apply(P[,16:18],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P$pat <- apply(P[,c(16,17,19)],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P <- P[,1:14]
write.table(P,"m116_chr11",sep="\t",row.names=F,quote=F)
#m116chr13
X <- read.table("m116_chr11.txt",header=T,sep="\t") 
P <- merge(X[X$phase=="PHASED",],Y[(Y$chr=="chr11"),],by.x="pos",by.y="pos")
P$mat <- apply(P[,16:18],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P$pat <- apply(P[,c(16,17,19)],1,function(x) if(x[3]=="1"){return(x[2])}else{return(x[1])})
P <- P[,1:14]
write.table(P,"m116_chr11",sep="\t",row.names=F,quote=F)
