### Sample plan data cleaning ###
library(stringr)
library(tidyverse)
library(dplyr)

#setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing/SamplingPlan")

#############################################
##### formatting the Sampling Plan data #####
#############################################
dat = read.csv(file = "../2016_Clayton/Field Book (2016) - Clayton - Sampling Plan_TIDIED.csv", header = TRUE)


#filter out the recongizable MO###
dat = dat %>% filter(str_detect(dat$Genotype, "Mo") == TRUE)
#create the breedtype category
BreedType = ifelse(substr(dat$Genotype, 1,1)=="M", "Inbred", "Hybrid")
dat = cbind(BreedType, dat)


snpFull = read.csv(file = "../IBM94markerset08seq.csv", header = TRUE)
snp = snpFull[,-(1:5)]

#############################################
### matching the Mo### to their snp Values ##
#############################################
library(data.table)
#function takes the last n characters of a string
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
#need = c(1,5,7,8,9,10,16) #column numbers of needed variables
#dat.num = dat.num[,need]
colnames(snpMatch) = "GenotypeNum"
#datnames =names(dat.num)
colnames(dat.num)[16] = "GenotypeNum"
dat2 = merge(dat.num,snpMatch, by.x = "GenotypeNum", by.y = "GenotypeNum",all.MoNum = all)
dat2 = dat2[order(dat2$Genotype, decreasing = FALSE),]
dat2 = dat2[,-1] #omit GenotypeNum
dat2[1:20, 1:20]
#write.csv(dat2, "SamplingPlan_dat2.csv")

########################################
######### Including plate data #########
########################################
colnames(dat2)[4] = "sampleID"

# 2018_11 plate data #
### the file in plate will come from the qPCR output including ratios ###
plate_11 = read.csv(file = "../2018_11/2018_11_withStress.csv", header = TRUE)
colnames(plate_11)[2] = "sampleID"
full_11 = merge(plate_11, dat2, by = "sampleID")

# 2018_6 plate data #
plate_6 = read.csv(file = , header = TRUE)
colnames(plate_11)[2] = "sampleID"
full_6 = merge(plate_11, dat2, by = "sampleID")

# 2018_8 plate data #
plate_8 = read.csv(file = , header = TRUE)
colnames(plate_11)[2] = "sampleID"
full_8 = merge(plate_11, dat2, by = "sampleID")

#####Adding marker location and chromosome#####

addmarker <- function(full, plate){
  aux = matrix(snpFull$incre_new, nrow= 1)
  aux = rbind(aux,snpFull$Chromosome)
  zeros = dim(full)[2]-dim(aux)[2]
  fillnames = names(full)[1:zeros]
  other = as.data.frame(matrix(rep(0,2*zeros), nrow = 2)) #repeat the number of 0 as the number of variables
  aux = cbind(other,aux)
  colnames(aux) = rep("",dim(aux)[2])
  colnames(full) = rep("",dim(aux)[2])
  dat3 = rbind(aux,full)
  colnames(dat3) = c(fillnames, as.character(snpFull$markername))
  return(dat3)
  
}

vqtl_11 = addmarker(full_11)
#write.csv(vqtl_11, file = "../2018_11/vqtlinput_11.csv" ,row.names = FALSE)
vqtl_6 = addmarker(full_6)
#write.csv(vqtl_11, file = "../2018_6/vqtlinput_6.csv" ,row.names = FALSE)
vqtl_8 = addmarker(full_8)
#write.csv(vqtl_11, file = "../2018_8/vqtlinput_8.csv" ,row.names = FALSE)

#####MAKE SURE TO DELETE THE EXTRA ZEROS IN [1:2,1:4] IN EXCEL AFTERWARDS#####

