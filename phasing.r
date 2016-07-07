options(stringsAsFactors = FALSE)
snps <- read.table("snps",header=T,sep="\t")
 
m27_trio <- snps[,c(1:6,13:14,7:8)]

m27_trio <- remove_impossible(m27_trio)
m27_trio <- remove_hom_child(m27_trio)
m27_trio <- remove_all_het(m27_trio)
phased_m27 <- cbind(m27_trio[,1:4],t(as.data.frame(phaser(m27_trio[,5:10]))))
colnames(phased_m27) <- c("chr","pos","ref","alt","m9.m27","m13.m27")

m116_trio <- snps[,c(1:4,11:12,7:8,9:10)]
m116_trio <- remove_impossible(m116_trio)
m116_trio <- remove_hom_child(m116_trio)
m116_trio <- remove_all_het(m116_trio)
phased_m116 <- cbind(m116_trio[,1:4],t(as.data.frame(phaser(m116_trio[,5:10]))))
colnames(phased_m116) <- c("chr","pos","ref","alt","mm106.m116","m27.m116")

all_phased <- merge(phased_m116,phased_m27,all=T)
all_phased[is.na(all_phased)] <- 0
# chr,pos,ref,alt,mm106.m116,m27.m116,m9.m27,m13.m27

##create VCF output
to_vcf <- all_phased
to_vcf$ID <- "."
to_vcf$QUAL <- "."
to_vcf$FILTER <- "PASS"
to_vcf$INFO <- "AF=1"
to_vcf$FORMAT <- "GT"
# chr,pos,ref,alt,mm106.m116,mm27.m116,m9.m27,m9.m27
to_vcf$m116 <- paste(to_vcf$mm106.m116,"|",to_vcf$m27.m116,sep="")
to_vcf$m27 <- paste(to_vcf$m9.m27,"|",to_vcf$m9.m27,sep="")
vcf2 <- to_vcf[,c(1,2,9,3,4,10:15)]
write('##fileformat=VCFv4.1',file="m116_m27.vcf",append=F)
write('##filedate=20160526',file="m116_m27.vcf",append=T)
write('##source="beagle.jar (r1399)"',file="m116_m27.vcf",append=T)
write('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated Allele Frequencies">',file="m116_m27.vcf",append=T)
write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',file="m116_m27.vcf",append=T)
write(paste("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","m116","m27",sep="\t"),file="m116_m27.vcf",append=T)
write.table(vcf2,file="m116_m27.vcf",sep="\t",quote=F,col.names=F,row.names=F,append=T)
###end vcf

# phase m116 utilising phased m27 data
m116_trio <- snps[,c(1:4,11:12,9:10)]
m116_trio <- merge(phased_m27,m116_trio,all.x=T)
m116_trio <- remove_impossible(m116_trio)
m116_trio <- remove_hom_child(m116_trio)
m116_trio <- remove_all_het(m116_trio)
phased_m116 <- t(as.data.frame(phaser(m116_trio[,5:10])))
fullphase<- cbind(m116_trio[,1:6],phased_m116[,1])
colnames(fullphase)[7] <- "m116"
 
remove_impossible <- function(X) {
	X <- X[
		(
			(X[,9]==X[,5]| X[,9]==X[,6])|
			(X[,10]==X[,5]| X[,10]==X[,6])
		) &
		( 	
			(X[,9]==X[,7]| X[,9]==X[,8]) |
			(X[,10]==X[,7]| X[,10]==X[,8])
		)
	,]
	return(X)
}

remove_hom_child <- function(X) {
	X <- X[
		X[,9]!=X[,10]
	,]
	return(X)
}

remove_all_het <- function(X) {
	X <- X[
		(X[,5]==X[,6])|
		(X[,7]==X[,8])|
		(X[,9]==X[,10])
	,]
	return(X)
}

phaser <- function(X) {
	apply(X,1,function(x) 
	{
		count=0
		if(x[5]==x[4]){count=count+1}
		if(x[5]==x[3]){count=count+2}
		if(x[5]==x[2]){count=count+4}
		if(x[5]==x[1]){count=count+8}
		if(x[5]==1) {
			if(count==1|count==2|count==3|count==7|count==11){
				return(c(0,1))
			} else {
				return(c(1,0))
			}
		} else {
			if(count==1|count==2|count==3|count==7|count==11){
				return(c(1,0))
			} else {
				return(c(0,1))
			}
		}
	})
}







