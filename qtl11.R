### four plots 
### m9 x m13
### m27
### mm106
### m116

X <- DW2

qtl11_markers <- read.table("qtl11_markers",header=T,sep="\t")
temp <- merge(qtl11_markers,snps[snps$chr=="chr11",],by.x="pos", by.y="pos")
temp$m9 <- rowSums(temp[,6:7])
temp$m13 <- rowSums(temp[,14:15])
temp$mm106 <- rowSums(temp[,10:11])

qtl11_usable <- read.table("q11_usable",header=T,sep="\t")


### m9 x m13
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
#axis(side=2,pos=1.1,lwd=2,at=ticks)

#m9
segments(1.7,min(X[,2]),1.7,max(X[,2]),col="gray75",lwd=6)
apply(temp[temp$m9>=1,],1,function(x) segments(1.725,as.numeric(x[1]),1.675,as.numeric(x[1]),col=1,lwd=1))
apply(qtl11_usable[qtl11_usable$m27=="m9",],1,function(x) segments(1.725,as.numeric(x[1]),1.675,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

#m13
segments(2.0,min(X[,2]),2.0,max(X[,2]),col="gray75",lwd=6)
apply(temp[temp$m13>=1,],1,function(x) segments(1.975,as.numeric(x[1]),2.025,as.numeric(x[1]),col=1,lwd=1))
apply(qtl11_usable[qtl11_usable$m27=="m13",],1,function(x) segments(1.975,as.numeric(x[1]),2.025,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

dev.off()

### m27
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")

#m9
segments(1.8,min(X[,2]),1.8,max(X[,2]),col="gray75",lwd=3)
apply(qtl11_usable[qtl11_usable$m27=="m9",],1,function(x) segments(1.81,as.numeric(x[1]),1.79,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

#m13
segments(2.0,min(X[,2]),2.0,max(X[,2]),col="gray75",lwd=3)
apply(qtl11_usable[qtl11_usable$m27=="m13",],1,function(x) segments(1.99,as.numeric(x[1]),2.01,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

dev.off()

### mm106

plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
segments(1.7,min(X[,2]),1.7,max(X[,2]),col="gray75",lwd=6)
apply(temp[temp$mm106>=1,],1,function(x) segments(1.725,as.numeric(x[1]),1.675,as.numeric(x[1]),col=1,lwd=0.5))

dev.off()


### m116
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")

segments(1.5,min(X[,2]),1.5,max(X[,2]),col="gray75",lwd=3)
apply(df11[df11$inht==2,],1,function(x) segments(1.8,as.numeric(x[1]),1.8,as.numeric(x[2]),col="gray75",lwd=3))
apply(df11[df11$inht==1,],1,function(x) segments(2,as.numeric(x[1]),2,as.numeric(x[2]),col="gray75",lwd=3))

apply(qtl11_usable[qtl11_usable$m27=="m13"&qtl11_usable$orig==1,],1,function(x) segments(1.99,as.numeric(x[1]),2.01,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))
apply(qtl11_usable[qtl11_usable$m27=="m9"&qtl11_usable$orig==2,],1,function(x) segments(1.81,as.numeric(x[1]),1.79,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

apply(qtl11_usable[qtl11_usable$m116=="mm106",],1,function(x) segments(1.485,as.numeric(x[1]),1.515,as.numeric(x[1]),col=1,lwd=0.5))

dev.off()


#segments(1.8,df11[22,1],1.8,df11[22,2],col="brown2",lwd=4,lend=2)
#segments(1.8,df11[26,1],1.8,df11[26,2],col="darkorchid1",lwd=4,lend=2)
#segments(1.8,df11[30,1],1.8,df11[30,2],col="darkseagreen1",lwd=4,lend=2)
#segments(1.8,df11[34,1],1.8,df11[34,2],col="darkslategray1",lwd=4,lend=2)
#segments(1.8,df11[46,1],1.8,df11[46,2],col="gold",lwd=4,lend=2)
