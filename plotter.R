### 1. this finds the number of matches between m116 c1 and m27.m9 in each run of 10 snps

### 2. then should look for stretches above a score of 8, then switch to m27.m13 when score goes below 2 
### and switch back when above 8 again

#1
phased_chr <- read.table("m116_phased_chr5.txt",header=T,sep="\t")
X <- phased_chr
X$tr <- (X[,4]==X[,3])
#cumulitive sum of differences between m116 c1 and m27.m9
csvec <- cumsum(X$tr)
# filter =  1  0  0  0  0  0  0  0  0  0 -1
f11 <- c(1,rep(0,9),-1)
# apply filter 
fv <- filter(csvec, f11, sides = 1)

#2
lapply(fv,function(x) if(x>=8){return(1)}else if(x<=2){return(2)}else{return(0)})
