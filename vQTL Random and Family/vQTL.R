install.packages("vqtl")
install.packages("qtl")
library("vqtl")
library("qtl")
family <-read.cross(dir = "~/PCR\ data/", file = "Family.csv")
family <- drop.nullmarkers(family)
#scan with variance
outv <-scanonevar(family)
family= sim.geno(family)
#trying to script fitqtls for all the markers
varlist<- paste("c",1:10, sep = "")
values <- rep(1,10)
df<-data.frame(varlist,values)
pval <- rep(2,3235)
names <- rep("empty",3235)
names[1:531] <- colnames(family$geno$`1`$data)
names[532:879] <- colnames(family$geno$`2`$data)
names[880:1263] <- colnames(family$geno$`3`$data)
names[1264:1583] <- colnames(family$geno$`4`$data)
names[1584:1908] <- colnames(family$geno$`5`$data)
names[1909:2205] <- colnames(family$geno$`6`$data)
names[2206:2475] <- colnames(family$geno$`7`$data)
names[2476:2734] <- colnames(family$geno$`8`$data)
names[2735:3003] <- colnames(family$geno$`9`$data)
names[3004:3235] <- colnames(family$geno$`10`$data)
for (x in 1:531) {
  for (y in 1:10) {
    if (length(family$geno[[y]]$data[1,]) >= x) {
      df[y,2] <- x
    }
    else {
      df[y,2] <- 1
    }
  }
  poslist<- c(outv[colnames(family$geno$`1`$data)[df[1,2]],2],
              outv[colnames(family$geno$`2`$data)[df[2,2]],2],
              outv[colnames(family$geno$`3`$data)[df[3,2]],2],
              outv[colnames(family$geno$`4`$data)[df[4,2]],2],
              outv[colnames(family$geno$`5`$data)[df[5,2]],2],
              outv[colnames(family$geno$`6`$data)[df[6,2]],2],
              outv[colnames(family$geno$`7`$data)[df[7,2]],2],
              outv[colnames(family$geno$`8`$data)[df[8,2]],2],
              outv[colnames(family$geno$`9`$data)[df[9,2]],2],
              outv[colnames(family$geno$`10`$data)[df[10,2]],2])
  qtl <- makeqtl(family, chr = 1:10, pos = poslist)
  fit <- fitqtl(family, qtl = qtl)
  sum <- summary(fit)
  pval[x] <- sum$result.drop[1,7]
  pval[df[2,2]+531] <- sum$result.drop[2,7]
  pval[df[3,2]+879] <- sum$result.drop[3,7]
  pval[df[4,2]+1263] <- sum$result.drop[4,7]
  pval[df[5,2]+1583] <- sum$result.drop[5,7]
  pval[df[6,2]+1908] <- sum$result.drop[6,7]
  pval[df[7,2]+2205] <- sum$result.drop[7,7]
  pval[df[8,2]+2475] <- sum$result.drop[8,7]
  pval[df[9,2]+2734] <- sum$result.drop[9,7]
  pval[df[10,2]+3003] <- sum$result.drop[10,7]
}
allpos<-outv[,2]
pval <- data.frame(names,allpos,pval)
pval
write.csv(pval, file = "FamilyPvalues.csv")

#Creating a Manhattan plot
#   Specifically using the qqman package
library(qqman)
CHR = rep(1,length(pval[,1]))
plotframe <- data.frame(SNP = pval[,1], BP = pval[,2], P = pval[,3], CHR = CHR)
manhattan(plotframe)
