### Combining Genomes and Mo### 

#setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing")
dat = read.csv(file = "Plant_Height.csv")
dat = dat[,-(4:5)]
#####Need to add in SNP info#####
snp = read.csv(file = "IBM94markerset08seq.csv")
snp = snp[,-(1:5)]
#snpT = t(snp)



relevant = data.frame(matrix(rep(0,length(dat$Genotype)*dim(snp)[1]), ncol = dim(snp)[1]))
dat2 = sapply(dat$Genotype, function(x){
  substr(x,2,2) = "O"
  column = which(colnames(snp) == x)
  vect = data.frame(as.character(snp[,column]))
  return(vect)
})
dat2 = as.data.frame(matrix(unlist(dat2), nrow = dim(dat)[1], byrow = TRUE))
colnames(dat2) = colnames()
library(beepr)
beep()
dim(dat3);dim(snp)
#####Adding back in the Trait info#####
dat3 = cbind(dat2$Height,dat2[,1:3],dat3)
colnames(dat3) = c(colnames(dat2[1:3]),"Height",as.character(snp$markername))
dat3[1:10,1:10]
#####Adding marker location and chromosome#####
aux = matrix(snp$incre_new, nrow= 1)
aux = rbind(aux,snp$Chromosome)
other = as.data.frame(matrix(rep(0,8), nrow = 2))
aux = cbind(other,aux)
colnames(aux) = rep("",3239)
colnames(dat3) = rep("",3239)
dat4 = rbind(aux,dat3)
colnames(dat4) = c("Height", colnames(dat2[1:3]),as.character(snp$markername))
dat4[1:10,1:10]
write.csv(dat4, file = ,
          row.names = FALSE)
beep()
#####MAKE SURE TO DELETE THE EXTRA ZEROS IN [1:2,1:4] IN EXCEL AFTERWARDS#####