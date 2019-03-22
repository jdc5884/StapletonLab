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
snp.num = as.integer(substrRight(as.character(colnames(snp)), 3))

snp.new = data.frame(lapply(snp,as.character),stringsAsFactors=FALSE)
snpMatch = rbind(snp.num,snp.new)
snpMatch = transpose(snpMatch)
# creating matching Genotype columns to merge the data
colnames(snpMatch) = c("GenotypeNum")
datnames =names(dat.num)[-7]
colnames(dat.num) = c(datnames,"GenotypeNum")
dat2 = merge(dat.num,snpMatch, by.x = "GenotypeNum", by.y = "GenotypeNum",all.MoNum = all)
dat2 = dat2[order(dat2$Genotype, decreasing = FALSE),]
dat2 = dat2[,-1]
write.csv(dat2, "SamplingPlan_dat2.csv")



