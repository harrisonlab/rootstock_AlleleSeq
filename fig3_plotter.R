options(stringsAsFactors = FALSE)

phased_chr <- read.table("m116_m13.m9",header=T,sep="\t")


# RBsnp1233 - RBbinsnp0048
DW1 <- phased_chr[phased_chr$chr=="chr5"&phased_chr$pos>=2426406&phased_chr$pos<=14674789,]
plotter_func(DW1,"chr5","DW1.pdf")
#  Roughly RBsnp1496 - Rbsnp1507 
DW2 <- phased_chr[phased_chr$chr=="chr11"&phased_chr$pos>=3470000&phased_chr$pos<=9965260,]
plotter_func(DW2,"chr11","DW2.pdf")
# RBsnp1597 - RBbinsnp0511
DW3 <- phased_chr[phased_chr$chr=="chr13"&phased_chr$pos>=34055&phased_chr$pos<=2813296,]
plotter_func(DW3,"chr13","DW3.pdf")

# font tester
#plot(NULL, xlim=c(0,4),ylim=c(0,length(names(pdfFonts()))),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
#sapply(names(pdfFonts()),function(q) text(x=rep(2,length(names(pdfFonts()))), y=which(names(pdfFonts())==q),pos=4, family=q,offset=0.1, cex=0.4,col=1,labels="This is a test"))


X <-DW1
chr <- "chr5"

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


dw5_genes <- read.table("dw5_genes",header=T,sep="\t")
dw5_genes <- dw5_genes[order(dw5_genes$start),]
dw5_genes$median <- apply(dw5_genes[,3:4],1,median)
dw5_genes$label <- 0

spacing <- (max(dw5_genes[dw5_genes$colour<=3,6]) - min(dw5_genes[dw5_genes$colour<=3,6])) / (29)
dw5_genes$x[dw5_genes$colour<=3] <-  c(1:30)
dw5_genes$x[dw5_genes$colour<=3] <-  (dw5_genes$x[dw5_genes$colour<=3]-1) *spacing +min(dw5_genes$median[dw5_genes$colour<=3])

spacing <- (max(dw5_genes[dw5_genes$colour>3,6]) - min(dw5_genes[dw5_genes$colour>3,6])) / (31)
dw5_genes$x[dw5_genes$colour>3] <-  c(1:length(dw5_genes$x[dw5_genes$colour>3]))
dw5_genes$x[dw5_genes$colour>3] <-  (dw5_genes$x[dw5_genes$colour>3]-1) *spacing +min(dw5_genes$median[dw5_genes$colour>3])

#pdf("test.pdf",width=8,height=8)  
ticks <- c(2400000,6400000,10400000,14400000)

plot(NULL, xlim=c(0,4),ylim=c(2400000,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
axis(side=2,pos=1.1,lwd=2,at=ticks)
apply(dw5_genes[dw5_genes$colour==1,],1,function(x) segments(2,as.numeric(x[6]),2.15,as.numeric(x[8]),col=1,lwd=1))
apply(dw5_genes[dw5_genes$colour==2,],1,function(x) segments(2,as.numeric(x[6]),2.15,as.numeric(x[8]),col=1,lwd=1))
apply(dw5_genes[dw5_genes$colour==3,],1,function(x) segments(2,as.numeric(x[6]),2.15,as.numeric(x[8]),col=1,lwd=1))

points(x=rep(2.15,length(dw5_genes$median[dw5_genes$colour==1])),cex=0.6, y=dw5_genes$x[(dw5_genes$colour==1)],pch=0,col=1)
points(x=rep(2.15,length(dw5_genes$median[dw5_genes$colour==2])),cex=0.6, y=dw5_genes$x[(dw5_genes$colour==2)],pch=1,col=1)
points(x=rep(2.15,length(dw5_genes$median[dw5_genes$colour==3])),cex=0.6, y=dw5_genes$x[(dw5_genes$colour==3)],pch=2,col=1)

apply(dw5_genes[dw5_genes$colour==4,],1,function(x) segments(1.9,as.numeric(x[6]),1.75,as.numeric(x[8]),col=1,lwd=1))
apply(dw5_genes[dw5_genes$colour==5,],1,function(x) segments(1.9,as.numeric(x[6]),1.75,as.numeric(x[8]),col=1,lwd=1))
apply(dw5_genes[dw5_genes$colour==6,],1,function(x) segments(1.9,as.numeric(x[6]),1.75,as.numeric(x[8]),col=1,lwd=1))

points(x=rep(1.75,length(dw5_genes$median[dw5_genes$colour==4])),cex=0.6, y=dw5_genes$x[(dw5_genes$colour==4)],pch=5,col=1)
points(x=rep(1.75,length(dw5_genes$median[dw5_genes$colour==5])),cex=0.6, y=dw5_genes$x[(dw5_genes$colour==5)],pch=6,col=1)
points(x=rep(1.75,length(dw5_genes$median[dw5_genes$colour==6])),cex=0.6, y=dw5_genes$x[(dw5_genes$colour==6)],pch=7,col=1)

apply(df[df$inht==1,],1,function(x) segments(1.9,as.numeric(x[1]),1.9,as.numeric(x[2]),col=1,lwd=3))
apply(df[df$inht==2,],1,function(x) segments(2,as.numeric(x[1]),2,as.numeric(x[2]),col=2,lwd=3))

text(x=rep(2.15,length(dw5_genes$median[dw5_genes$colour==1])), y=dw5_genes$x[(dw5_genes$colour==1)],family="Helvetica-Narrow",pos=4, offset=0.4, cex=0.6,col=1,labels=dw5_genes$Gene_ID[(dw5_genes$colour==1)])
text(x=rep(2.15,length(dw5_genes$median[dw5_genes$colour==2])), y=dw5_genes$x[(dw5_genes$colour==2)], family="Helvetica-Narrow",pos=4, offset=0.4, cex=0.6,col=1,labels=dw5_genes$Gene_ID[(dw5_genes$colour==2)])
text(x=rep(2.15,length(dw5_genes$median[dw5_genes$colour==3])), y=dw5_genes$x[(dw5_genes$colour==3)],family="Helvetica-Narrow",pos=4, offset=0.4, cex=0.6,col=1,labels=dw5_genes$Gene_ID[(dw5_genes$colour==3)])

text(x=rep(1.75,length(dw5_genes$median[dw5_genes$colour==4])), y=dw5_genes$x[(dw5_genes$colour==4)], family="Helvetica-Narrow",pos=2, offset=0.4, cex=0.6,col=1,labels=dw5_genes$Gene_ID[(dw5_genes$colour==4)])
text(x=rep(1.75,length(dw5_genes$median[dw5_genes$colour==5])), y=dw5_genes$x[(dw5_genes$colour==5)], family="Helvetica-Narrow",pos=2, offset=0.4, cex=0.6,col=1,labels=dw5_genes$Gene_ID[(dw5_genes$colour==5)])
text(x=rep(1.75,length(dw5_genes$median[dw5_genes$colour==6])), y=dw5_genes$x[(dw5_genes$colour==6)], family="Helvetica-Narrow",pos=2, offset=0.4, cex=0.6,col=1,labels=dw5_genes$Gene_ID[(dw5_genes$colour==6)])

dev.off()



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


dw11_genes <- read.table("dw11_genes",header=T,sep="\t")
dw11_genes <- dw11_genes[order(dw11_genes$start),]
dw11_genes$median <- apply(dw11_genes[,3:4],1,median)
dw11_genes$label <- 0

spacing <- (max(dw11_genes[dw11_genes$colour!=3,6]) - min(dw11_genes[dw11_genes$colour!=3,6])) / (41)
dw11_genes$x[dw11_genes$colour!=3] <-  c(1:42)
dw11_genes$x[dw11_genes$colour!=3] <-  (dw11_genes$x[dw11_genes$colour!=3]-1) *spacing +min(dw11_genes$median[dw11_genes$colour!=3])
dw11_genes$x[dw11_genes$colour==3] <- c( 3680700, 3850631, 9686142)

#pdf("test.pdf",width=8,height=8)  
ticks <- c(3000000,6000000,900000)

plot(NULL, xlim=c(0,4),ylim=c(3000000,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
axis(side=2,pos=1.1,lwd=2)

apply(dw11_genes[dw11_genes$colour==1,],1,function(x) segments(2,as.numeric(x[6]),2.15,as.numeric(x[8]),col=1,lwd=1))
apply(dw11_genes[dw11_genes$colour==2,],1,function(x) segments(2,as.numeric(x[6]),2.15,as.numeric(x[8]),col=1,lwd=1))
apply(dw11_genes[dw11_genes$colour==3,],1,function(x) segments(1.9,as.numeric(x[6]),1.75,as.numeric(x[8]),col=1,lwd=1))

points(x=rep(2.15,length(dw11_genes$median[dw11_genes$colour==1])),cex=0.6, y=dw11_genes$x[(dw11_genes$colour==1)],pch=15,col=1)
points(x=rep(2.15,length(dw11_genes$median[dw11_genes$colour==2])), cex=0.6,y=dw11_genes$x[(dw11_genes$colour==2)],pch=16,col=1)
points(x=rep(1.75,length(dw11_genes$median[dw11_genes$colour==3])),cex=0.6, y=dw11_genes$x[(dw11_genes$colour==3)],pch=17,col=1)

apply(df[df$inht==1,],1,function(x) segments(1.9,as.numeric(x[1]),1.9,as.numeric(x[2]),col=1,lwd=3))
apply(df[df$inht==2,],1,function(x) segments(2,as.numeric(x[1]),2,as.numeric(x[2]),col=2,lwd=3))

text(x=rep(2.15,length(dw11_genes$median[dw11_genes$colour==1])), y=dw11_genes$x[(dw11_genes$colour==1)],family="Helvetica-Narrow",pos=4, offset=0.4, cex=0.6,col=1,labels=dw11_genes$Gene_ID[(dw11_genes$colour==1)])
text(x=rep(2.15,length(dw11_genes$median[dw11_genes$colour==2])), y=dw11_genes$x[(dw11_genes$colour==2)],family="Helvetica-Narrow", pos=4, offset=0.4, cex=0.6,col=1,labels=dw11_genes$Gene_ID[(dw11_genes$colour==2)])
text(x=rep(1.75,length(dw11_genes$median[dw11_genes$colour==3])), y=dw11_genes$x[(dw11_genes$colour==3)],family="Helvetica-Narrow", pos=2, offset=0.4, cex=0.6,col=1,labels=dw11_genes$Gene_ID[(dw11_genes$colour==3)])

dev.off()

X <-DW3
chr <- "chr13"
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

dw13_genes <- read.table("dw13_genes",header=T,sep="\t")
dw13_genes <- dw13_genes[order(dw13_genes$start),]
dw13_genes$median <- apply(dw13_genes[,3:4],1,median)
dw13_genes$median <- floor(dw13_genes$median)
dw13_genes$label <- 0

dw13_genes[,6] <- c(245466,900578,1074849,1199142,2193540)

plot(NULL, xlim=c(0,4),ylim=c(0,max(X$pos)),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
axis(side=2,pos=1.1,lwd=2)

apply(dw13_genes,1,function(x) segments(1.9,as.numeric(x[6]),1.75,as.numeric(x[6]),col=1,lwd=1))

points(x=rep(1.75,length(dw13_genes[,6])),cex=1, y=dw13_genes[,6],pch=18,col=1)

apply(df[df$inht==1,],1,function(x) segments(1.9,as.numeric(x[1]),1.9,as.numeric(x[2]),col=1,lwd=3))
apply(df[df$inht==2,],1,function(x) segments(2,as.numeric(x[1]),2,as.numeric(x[2]),col=2,lwd=3))

text(x=rep(1.75,length(dw13_genes[,7])), y=dw13_genes[,6], family="Helvetica-Narrow",pos=2, offset=0.4, cex=0.6,col=1,labels=dw13_genes$Gene_ID)

dev.off()


