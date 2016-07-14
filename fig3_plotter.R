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

#points(rep(3.2,length(dw11_genes$median[(dw11_genes$colour==1)])),dw11_genes$median[(dw11_genes$colour==1)],col=5)
#labypos <- dw11_genes$median[(dw11_genes$colour==1)] + dw11_genes$label[dw11_genes$colour==1]


text(x=rep(2.4,length(dw11_genes$median[dw11_genes$colour==1])), y=dw11_genes$x[(dw11_genes$colour==1)],pos=4, offset=0.1, cex=0.4,col=4,labels=dw11_genes$Gene_ID[(dw11_genes$colour==1)])
text(x=rep(2.4,length(dw11_genes$median[dw11_genes$colour==2])), y=dw11_genes$x[(dw11_genes$colour==2)], pos=4, offset=0.1, cex=0.4,col=6,labels=dw11_genes$Gene_ID[(dw11_genes$colour==2)])

dev.off()
