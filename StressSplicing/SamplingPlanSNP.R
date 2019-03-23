### Sample plan data cleaning ###
library(stringr)
library(tidyverse)
library(dplyr)

#setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing")

dat = read.csv(file = "SamplingPlan.csv", header = TRUE)

dat = dat %>% filter(str_detect(dat$Genotype, "Mo") == TRUE)

BreedType = ifelse(substr(dat$Genotype, 1,1)=="M", "Inbred", "Hybrid")
dat = cbind(BreedType, dat)


snpFull = read.csv(file = "IBM94markerset08seq.csv", header = TRUE)
snp = snpFull[,-(1:5)]


library(data.table)
#function takes the last three values of a string
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
#creating datasets exclusively containing the last three numbers 
dat.num =cbind(dat, as.integer(substrRight(as.character(dat$Genotype), 3)))
#produces warning sign about NA's for data that is not availible for invalid Mo###
snp.num = as.integer(substrRight(as.character(colnames(snp)), 3))

snp.new = data.frame(lapply(snp,as.character),stringsAsFactors=FALSE)
snpMatch = rbind(snp.num,snp.new)
snpMatch = transpose(snpMatch)
# creating matching Genotype columns to merge the data
need = c(1,5,7,8,15) #column numbers of needed variables
dat.num = dat.num[,need]
colnames(snpMatch) = c("GenotypeNum")
datnames =names(dat.num)
colnames(dat.num) = c(datnames[-5],"GenotypeNum")
dat2 = merge(dat.num,snpMatch, by.x = "GenotypeNum", by.y = "GenotypeNum",all.MoNum = all)
dat2 = dat2[order(dat2$Genotype, decreasing = FALSE),]
dat2 = dat2[,-1]
dat2[1:20, 1:20]
write.csv(dat2, "SamplingPlan_dat2.csv")


#####Adding marker location and chromosome#####
aux = matrix(snpFull$incre_new, nrow= 1)
aux = rbind(aux,snpFull$Chromosome)
other = as.data.frame(matrix(rep(0,8), nrow = 2))
aux = cbind(other,aux)
colnames(aux) = rep("",dim(dat2)[2])
colnames(dat2) = rep("",dim(dat2)[2])
dat3 = rbind(aux,dat2)
colnames(dat3) = c("BreedType", "Genotype", "Derivation", "Barcode", as.character(snpFull$markername))
dat3[1:10,1:10]
write.csv(dat3, file = "snpSamplingPlan.csv" ,
          row.names = FALSE)
#####MAKE SURE TO DELETE THE EXTRA ZEROS IN [1:2,1:4] IN EXCEL AFTERWARDS#####


