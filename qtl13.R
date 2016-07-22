### four plots 
### m9 x m13
### m27
### mm106
### m116
X<-DW3
qtl13_markers <- read.table("qtl13_markers",header=T,sep="\t")
temp <- merge(qtl13_markers,snps[snps$chr=="chr13",],by.x="pos", by.y="pos")
temp$m9 <- rowSums(temp[,6:7])
temp$m13 <- rowSums(temp[,14:15])
temp$mm106 <- rowSums(temp[,10:11])


qtl13_usable <- read.table("q13_usable",header=T,sep="\t")


### m9 x m13
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
#axis(side=2,pos=1.1,lwd=2,at=ticks)

#m9
segments(1.7,min(X[,2]),1.7,max(X[,2]),col="gray75",lwd=6)
apply(temp[temp$m9>=1,],1,function(x) segments(1.725,as.numeric(x[1]),1.675,as.numeric(x[1]),col=1,lwd=0.5))
apply(qtl13_usable[qtl13_usable$m27=="m9",],1,function(x) segments(1.725,as.numeric(x[1]),1.675,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

#m13
segments(2.0,min(X[,2]),2.0,max(X[,2]),col="gray75",lwd=6)
apply(temp[temp$m13>=1,],1,function(x) segments(1.975,as.numeric(x[1]),2.025,as.numeric(x[1]),col=1,lwd=0.5))
apply(qtl13_usable[qtl13_usable$m27=="m13",],1,function(x) segments(1.975,as.numeric(x[1]),2.025,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

dev.off()

### m27
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")

#m9
segments(1.8,min(X[,2]),1.8,max(X[,2]),col="gray75",lwd=3)
apply(qtl13_usable[qtl13_usable$m27=="m9",],1,function(x) segments(1.81,as.numeric(x[1]),1.79,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

#m13
segments(2.0,min(X[,2]),2.0,max(X[,2]),col="gray75",lwd=3)
apply(qtl13_usable[qtl13_usable$m27=="m13",],1,function(x) segments(1.99,as.numeric(x[1]),2.01,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

dev.off()

### mm106
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")

segments(1.7,min(X[,2]),1.7,max(X[,2]),col="gray75",lwd=6)
apply(temp[temp$mm106>=1,],1,function(x) segments(1.725,as.numeric(x[1]),1.675,as.numeric(x[1]),col=1,lwd=0.5))

dev.off()


### m116
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")

segments(1.5,min(X[,2]),1.5,max(X[,2]),col="gray75",lwd=3)
apply(df13[df13$inht==2,],1,function(x) segments(1.8,as.numeric(x[1]),1.8,as.numeric(x[2]),col="gray75",lwd=3))
apply(df13[df13$inht==1,],1,function(x) segments(2,as.numeric(x[1]),2,as.numeric(x[2]),col="gray75",lwd=3))

apply(qtl13_usable[qtl13_usable$m27=="m9"&qtl13_usable$orig==2,],1,function(x) segments(1.81,as.numeric(x[1]),1.79,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

apply(qtl13_usable[qtl13_usable$m27=="m13"&qtl13_usable$orig==1,],1,function(x) segments(1.99,as.numeric(x[1]),2.01,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))
apply(qtl13_usable[qtl13_usable$m116=="mm106"|qtl13_usable$m116=="11",],1,function(x) segments(1.485,as.numeric(x[1]),1.515,as.numeric(x[1]),col=1,lwd=0.5))

dev.off()


segments(1.8,df13[8,1],1.8,df13[8,2],col="coral1",lwd=3,lend=2)
#segments(1.8,df13[12,1],1.8,df13[12,2],col="darkgreen",lwd=3,lend=2)
#segments(1.8,df13[16,1],1.8,df13[16,2],col="blue3",lwd=3,lend=2)

#text(x=rep(2.15,length(qtl13_usable$m27=="m13")), y=qtl13_usable$x[qtl13_usable$m27=="m13"],family="Helvetica-Narrow",pos=4, offset=0.4, cex=0.6,col=1,labels=qtl13_usable$name[qtl13_usable$m27=="m13"])
#text(x=rep(1.35,length(qtl13_usable[qtl13_usable$m116=="mm106",])), y=qtl13_usable$x[qtl13_usable$m116=="mm106"],family="Helvetica-Narrow", pos=2, offset=0.4, cex=0.6,col=1,labels=qtl13_usable$name[qtl13_usable$m116=="mm106"])
