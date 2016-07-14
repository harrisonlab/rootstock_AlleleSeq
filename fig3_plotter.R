phased_chr <- read.table("m116_m13.m9",header=T,sep="\t")
DW2 <- phased_chr[phased_chr$chr=="chr11"&phased_chr$pos>=3470000&phased_chr$pos<=9965260,]
X <-DW2
chr <- "chr11"

myc <- X[X$chr==chr,]
# order by snp position
myc<-myc[order(myc$pos),]
# Test whether the snp in m116 was inherited from m13 (via m27) 
myc$tr <- (myc[,7]==myc[,6])
# cumulitive sum of differences between m116 c5 and m27.m13
csvec <- cumsum(myc$tr)
# filter =  1  0  0  0  0  0  0  0  0  0 -1
f11 <- c(1,rep(0,19),-1)
# apply filter 
fv <- filter(csvec, f11, sides = 1)
# first 9 entries must be NA, set them to 5 
fv[is.na(fv)] <-10

#2 get +15 sequences
mv <- 0
myc$new <- sapply(fv,function(x) if(x>=15){mv<<-1;return(mv)}else if(x<=5){mv<<-2;return(mv)}else if(x==10){mv<<-0;return(mv)}else{return(mv)})
# find postion where sequence changes
change <- myc[c(1,1+which(diff(myc$new)!=0)),]
# get stop and start positions from change - can't think of a good R way of doing this so will use a for loop
df <- data.frame(start=as.integer(),end=as.integer(),inht=as.integer())
for (i in 1:(length(change[,1])-1)) {
  start <- change[i,2]
  end <- change[i+1,2]
  inht <- change[i,9]
  df[nrow(df)+1,] <- c(start,end,inht)
}

dw11_genes <- read.table("dw11_genes",header=T\sep="\t")
dw11_genes$median <- apply(dw11_genes[,3:4],1,median)
dw11_genes$label <- 0
spacing <- (max(dw11_genes[,6]) -min(dw11_genes[,6])) / (length(dw11_genes$median)-1)
dw11_genes$x <-  c(1:42)
dw11_genes$x <-  (dw11_genes$x-1) *spacing +min(dw11_genes$median)


plot(NULL, xlim=c(0,4),ylim=c(3000000,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
axis(side=2,pos=1.6,lwd=2)
apply(dw11_genes[dw11_genes$colour==1,],1,function(x) segments(2.2,as.numeric(x[6]),2.4,as.numeric(x[8]),col=4,lwd=1))
apply(dw11_genes[dw11_genes$colour==2,],1,function(x) segments(2.2,as.numeric(x[6]),2.4,as.numeric(x[8]),col=6,lwd=1))

apply(df[df$inht==1,],1,function(x) segments(1.9,as.numeric(x[1]),1.9,as.numeric(x[2]),col=1,lwd=3))
apply(df[df$inht==2,],1,function(x) segments(2.2,as.numeric(x[1]),2.2,as.numeric(x[2]),col=2,lwd=3))

text(x=rep(2.4,length(dw11_genes$median[dw11_genes$colour==1])), y=dw11_genes$x[(dw11_genes$colour==1)],pos=4, offset=0.1, cex=0.4,col=4,labels=dw11_genes$Gene_ID[(dw11_genes$colour==1)])
text(x=rep(2.4,length(dw11_genes$median[dw11_genes$colour==2])), y=dw11_genes$x[(dw11_genes$colour==2)], pos=4, offset=0.1, cex=0.4,col=6,labels=dw11_genes$Gene_ID[(dw11_genes$colour==2)])

dev.off()
