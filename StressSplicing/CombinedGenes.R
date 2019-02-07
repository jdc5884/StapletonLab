### Combining Genomes and Mo### 
library(stringr)
library(tidyverse)
library(dplyr)

#setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing")
dat = read.csv(file = "Plant_Height_Update.csv", header = TRUE)

#Take out unneeded IBMB###, NA, B73 loci
dat = dat %>% filter(str_detect(dat$Height..in.., "N") == FALSE)
dat = dat %>% filter(str_detect(dat$Genotype, "B") == FALSE)
dat$Genotype = str_remove(dat$Genotype, " ") #removes any empty space in Mo###

#Create Categorical Variables for PH207*Mo### and Mo### by gene breed
BreedType = ifelse(substr(dat$Genotype, 1,1)=="M", "Inbred", "Outbred")
dat = cbind(dat, BreedType)

#Add in SNP info from Marker data CSV, beginning with column six
snpFull = read.csv(file = "IBM94markerset08seq.csv", header = TRUE)
snp = snpFull[,-(1:5)]

### Assigning Genotypes to Mo###
library(data.table)
#function takes the last three values of a string
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
#creating datasets exclusively containing the last three numbers 
dat.num = substrRight(as.character(dat$Genotype), 3)
snp.num = substrRight(as.character(colnames(snp)), 3)

snp.new = data.frame(lapply(snp,as.character),stringsAsFactors=FALSE)
snpMatch = rbind(snp.num,snp.new)
snpMatch = transpose(snpMatch)
# creating matching Genotype columns to merge the data
colnames(snpMatch) = c("Genotype")
MoNum = cbind(dat.num)
colnames(MoNum) = c("Genotype")
dat2 = merge(MoNum,snpMatch,by = "Genotype")
dat2[,1:10]
dat2 = dat2[,-1]
#combine with the full data description
dat2 = cbind(dat, dat2)
write.csv(dat2, "dat2.csv")



## Issues here ##
#####################################################################################
#####Adding marker location and chromosome#####
aux = matrix(snpFull$incre_new, nrow= 1)
aux = rbind(aux,snpFull$Chromosome)
other = as.data.frame(matrix(rep(0,8), nrow = 2))
aux = cbind(other,aux)
colnames(aux) = rep("",3239)
colnames(dat2) = rep("",3239)
dat3 = rbind(aux,dat2) #dimentions do not match to bind both dataframes
colnames(dat3) = c("Height", colnames(dat2[1:3]),as.character(snp$markername))
dat3[1:10,1:10]
write.csv(dat3, file = ,
          row.names = FALSE)
beep()
#####MAKE SURE TO DELETE THE EXTRA ZEROS IN [1:2,1:4] IN EXCEL AFTERWARDS#####