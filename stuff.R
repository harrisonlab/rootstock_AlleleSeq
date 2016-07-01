# shell stuff
# grep ^5 03-RNA_L1_1.fq.trim.counts.txt >S3_chr5
# grep ^5 05-RNA_L1_1.fq.trim.counts.txt >S5_chr5
# grep ^5 07-RNA_L1_1.fq.trim.counts.txt >S7_chr5
options(stringsAsFactors = FALSE)

args <- commandArgs(TRUE) 
args[1] <- "S9_chr5"
args[2] <- "S11_chr5"
args[3] <- "S13_chr5"

t1 <- read.table(args[1],header=F,sep="\t")
t2 <- read.table(args[2],header=F,sep="\t")
t3 <- read.table(args[3],header=F,sep="\t")
x1 <- merge(t1[,c(2,7:16)],t2[,c(2,7:16)],all=T,by.x="V2",by.y="V2")
x2 <- merge(x1,t3[,c(2,7:16)],all=T,by.x="V2",by.y="V2")

myfunc <- function(v) {
  n=""
  if(!is.na(v[1])) {n=v[1]}
  else if(!is.na(v[2])) {n=v[2]}
  else {n=v[3]}
  return(n)
}
x2$pos <- x2[,1]
x2$phase <- apply(x2[,c(2,12,22)],1,function(x) myfunc(x))  
x2$mat <- apply(x2[,c(3,13,23)],1,function(x) myfunc(x))
x2$pat <- apply(x2[,c(4,14,24)],1,function(x) myfunc(x))

x2[is.na(x2)] <- 0

x2$A_sum <- rowSums(x2[,c(5,15,25)]) 
x2$C_sum <- rowSums(x2[,c(6,16,26)])
x2$G_sum <- rowSums(x2[,c(7,17,27)])
x2$T_sum <- rowSums(x2[,c(8,18,28)])


x2$S3_dtl <- paste(x2[,10],x2[,5],x2[,6],x2[,7],x2[,8],sep=":") # Sx..
x2$S5_dtl <- paste(x2[,20],x2[,15],x2[,16],x2[,17],x2[,18],sep=":")
x2$S7_dtl <- paste(x2[,30],x2[,25],x2[,26],x2[,27],x2[,28],sep=":")
x2$S3_padj <- x2[,11]
x2$S5_padj <- x2[,21]
x2$S7_padj <- x2[,31]

x3 <- x2[,32:45]
