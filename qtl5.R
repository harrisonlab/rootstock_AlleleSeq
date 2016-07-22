### four plots 
### m9 x m13
### m27
### mm106
### m116

X <- DW1

qtl5_markers <- read.table("qtl5_markers",header=T,sep="\t")
temp <- merge(qtl5_markers[,1:2],snps[snps$chr=="chr5",],by.x="pos", by.y="pos")
temp$m9 <- rowSums(temp[,6:7])
temp$m13 <- rowSums(temp[,14:15])
temp$mm106 <- rowSums(temp[,10:11])

qtl5_usable <- read.table("q5_usable",header=T,sep="\t")


### m9 x m13
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
#axis(side=2,pos=1.1,lwd=2,at=ticks)

#m9
segments(1.7,min(X[,2]),1.7,max(X[,2]),col="gray75",lwd=6)
apply(temp[temp$m9>=1,],1,function(x) segments(1.725,as.numeric(x[1]),1.675,as.numeric(x[1]),col=1,lwd=1))
apply(qtl5_usable[qtl5_usable$m27=="m9",],1,function(x) segments(1.725,as.numeric(x[1]),1.675,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

#m13
segments(2.0,min(X[,2]),2.0,max(X[,2]),col="gray75",lwd=6)
apply(temp[temp$m13>=1,],1,function(x) segments(1.975,as.numeric(x[1]),2.025,as.numeric(x[1]),col=1,lwd=1))
#apply(qtl5_usable[qtl5_usable$m27=="m13"&qtl5_usable$orig==1,],1,function(x) segments(1.99,as.numeric(x[1]),2.01,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))
apply(qtl5_usable[qtl5_usable$m27=="m13",],1,function(x) segments(1.975,as.numeric(x[1]),2.025,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

dev.off()


### m27
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")

#m9
segments(1.8,min(X[,2]),1.8,max(X[,2]),col="gray75",lwd=3)
apply(qtl5_usable[qtl5_usable$m27=="m9",],1,function(x) segments(1.81,as.numeric(x[1]),1.79,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

#m13
segments(2,min(X[,2]),2,max(X[,2]),col="gray75",lwd=3)
#apply(qtl5_usable[qtl5_usable$m27=="m13",],1,function(x) segments(1.985,as.numeric(x[1]),2.015,as.numeric(x[1]),col=1,lwd=0.5))
apply(qtl5_usable[qtl5_usable$m27=="m13",],1,function(x) segments(1.99,as.numeric(x[1]),2.01,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))

dev.off()

### mm106

plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")
segments(1.7,min(X[,2]),1.7,max(X[,2]),col="gray75",lwd=6)
apply(temp[temp$mm106>=1,],1,function(x) segments(1.725,as.numeric(x[1]),1.675,as.numeric(x[1]),col=1,lwd=0.5))

dev.off()


### m116
plot(NULL, xlim=c(0,4),ylim=c(0,max(X[,2])),xlab = "",tck=0,xaxt='n',yaxt='n',ann=F,bty="n")

segments(1.5,min(X[,2]),1.5,max(X[,2]),col="gray75",lwd=3)
apply(df5[df5$inht==2,],1,function(x) segments(1.8,as.numeric(x[1]),1.8,as.numeric(x[2]),col="gray75",lwd=3))
apply(df5[df5$inht==1,],1,function(x) segments(2,as.numeric(x[1]),2,as.numeric(x[2]),col="gray75",lwd=3))

apply(qtl5_usable[qtl5_usable$m27=="m13"&qtl5_usable$orig==1,],1,function(x) segments(1.99,as.numeric(x[1]),2.01,as.numeric(x[1]),col=x[6],lwd=1.5,lend=2))
apply(qtl5_usable[qtl5_usable$m116=="mm106",],1,function(x) segments(1.485,as.numeric(x[1]),1.515,as.numeric(x[1]),col=1,lwd=0.5))

dev.off()

#segments(2,14520884,2,14520985,col="pink",lwd=6,lend=2) #0047
#segments(2,8643407,2,8643407,col="orange",lwd=3,lend=2) # 0042
#segments(2,14365897,2,14438057,col="orange",lwd=3,lend=2) # 0042
#segments(2,11970050,2,11970050,col="green",lwd=3,lend=2) # 0045
#segments(2,11985564,2,12016362,col="green",lwd=3,lend=2) # 0045
#segments(2,11479429,2,11479429,col="cyan",lwd=3,lend=2) # 0044
#segments(2,9491908,2,9512695,col="mediumblue",lwd=3,lend=2) # 0043

#segments(2,df5[62,1],2,df5[62,2],col="pink",lwd=3,lend=2)
#segments(2,df5[46,1],2,df5[46,2],col="orange",lwd=3,lend=2)
#segments(2,df5[42,1],2,df5[42,2],col="green",lwd=3,lend=2)
#segments(2,df5[40,1],2,df5[40,2],col="cyan",lwd=3,lend=2)