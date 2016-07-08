### 1. this finds the number of matches between m116 c1 and m27.m9 in each run of 20 snps

### 2. then should look for stretches above a score of 8 (i.e. inheritted from m13, then switch to m27.m9 when score goes below 2 
### and switch back when above 8 again

### 3. plot each block of inherited snps using plot segments

#1
phased_chr <- read.table("m116_m13.m9",header=T,sep="\t")

plotter_func <- function(X,chr) {
	# get chromosome
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

	mv <- 0
	#2 get +8 sequences
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

	#3 plotting
	pdf(paste(chr,".pdf",sep=""),width=1.5,height=6)  
	plot(NULL, xlim=c(0.8,1.3),ylim=c(0,max(X[,2])),xlab = "",tck=0,yaxt='n',ann=F,bty="n")
	#,at=c(1,2,3,4),labels=c("ca1_p1","ca1_p2","ca2_p1","ca2_p2")
	apply(df[df$inht==1,],1,function(x) segments(0.9,as.numeric(x[1]),0.9,as.numeric(x[2]),col=1,lwd=3))
	apply(df[df$inht==2,],1,function(x) segments(1.2,as.numeric(x[1]),1.2,as.numeric(x[2]),col=2,lwd=3))
	dev.off()
}

for (i in 1:17) {
	plotter_func(phased_chr,paste("chr",i,sep=""))
}
