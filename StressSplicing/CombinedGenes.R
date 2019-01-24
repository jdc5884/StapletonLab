### Combining Genomes and Mo### 
library(stringr)
library(tidyverse)

#setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing")
dat = read.csv(file = "Plant_Height.csv", header = TRUE)

#Take out unneeded IBMB###, NA, B73 loci
dat = dat %>% filter(str_detect(dat$Height..in.., "N") == FALSE)
dat = dat %>% filter(str_detect(dat$Genotype, "B") == FALSE)

#Create Categorical Variables for PH207*Mo### and Mo### by gene breed
BreedType = ifelse(substr(dat$Genotype, 1,1)=="M", "Inbred", "Outbred")
dat = cbind(dat, BreedType)

#Add in SNP info from Marker data CSV, beginning with column six
snp = read.csv(file = "IBM94markerset08seq.csv", header = TRUE)
snp = snp[,-(1:5)]

### Assigning Genotypes to Mo###
### Austin's code ###
library(data.table)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

dat.num = substrRight(as.character(dat$Genotype), 3) # Leaves out begining 0 for Mo010 and Mo039 
snp.num = substrRight(as.character(colnames(snp)), 3)

snp.new = data.frame(lapply(snp,as.character),stringsAsFactors=FALSE)
snpMatch = rbind(snp.num,snp.new)
snpMatch = transpose(snpMatch)
# creating matching Genotype columns to merge the data
colnames(snpMatch) = c("Genotype")
MoNum = cbind(dat.num)
colnames(MoNum) = c("Genotype")
dat2 = merge(MoNum,snpMatch,by = "Genotype") ### issue: the dim of dat2 should = MoNum when it is merged with the Genotypes (1527-1437)
head(dat2[,1:10],20)
# this code puts in extra Mo### values when merging the two sets

### end Austin's code ###

#Running into errors beginning with colnames(dat2)#
colnames = colnames(dat2)
library(beepr)
beep()
dim(dat2);dim(snp)

#####Adding back in the Trait info#####
dat2 = cbind(dat$Height,dat[,1:4],dat2)
colnames(dat2) = c(colnames(dat[1:4]),"Height",as.character(snp$markername))
dat2[1:10,1:10]




## Issues here ##
#####################################################################################
#####Adding marker location and chromosome#####
aux = matrix(snp$incre_new, nrow= 1)
aux = rbind(aux,snp$Chromosome)
other = as.data.frame(matrix(rep(0,8), nrow = 2))
aux = cbind(other,aux)
colnames(aux) = rep("",3239)
colnames(dat2) = rep("",3239)
dat3 = rbind(aux,dat2)
colnames(dat3) = c("Height", colnames(dat2[1:3]),as.character(snp$markername))
dat3[1:10,1:10]
write.csv(dat3, file = ,
          row.names = FALSE)
beep()
#####MAKE SURE TO DELETE THE EXTRA ZEROS IN [1:2,1:4] IN EXCEL AFTERWARDS#####