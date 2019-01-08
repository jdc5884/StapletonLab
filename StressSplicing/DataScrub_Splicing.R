#####Scrubbing Manching Combined Stress Data#####
dat = read.csv(file = "C:/Users/twili/Desktop/GIThub/Stapleton-Lab/Manching BayesNet/Manching2012PlantHT.csv")
unique(dat$Line)
#Standardizing the 11 unique treatment descriptions to the correct number 8
unique(dat$Treatment)
unique(dat$Treatment)[2] == unique(dat$Treatment)[3]
levels(dat$Treatment)[3] = levels(dat$Treatment)[2]
levels(dat$Treatment)[10] = levels(dat$Treatment)[9]
levels(dat$Treatment)[4] = levels(dat$Treatment)[8]
length(levels(dat$Treatment)) == 8
levels(dat$Treatment) = c("control","ln","lw,ln","lw,p",
                          "ln,p","lw,ln,p","lw","p")
levels(dat$Treatment)
#Setting up quantitative columns using treatment levels
treats = matrix(rep(0,length(dat$Treatment) * 3), ncol = 3)
treats = t(sapply(1:length(dat$Treatment), function(x){
  if(as.character(dat$Treatment[x]) == "control"){
    treats[x,] = c(0,0,0)
  }else if(as.character(dat$Treatment[x]) == "ln"){
    treats[x,] = c(0,1,0)
  }else if(as.character(dat$Treatment[x]) == "lw,ln"){
    treats[x,] = c(1,1,0)
  }else if(as.character(dat$Treatment[x]) == "lw,p"){
    treats[x,] = c(1,0,1)
  }else if(as.character(dat$Treatment[x]) == "ln,p"){
    treats[x,] = c(0,1,1)
  }else if(as.character(dat$Treatment[x]) == "lw,ln,p"){
    treats[x,] = c(1,1,1)
  }else if(as.character(dat$Treatment[x]) == "lw"){
    treats[x,] = c(1,0,0)
  }else if(as.character(dat$Treatment[x]) == "p"){
    treats[x,] = c(0,0,1)}
}))
treats = as.data.frame(treats)
colnames(treats) = c("Low Water", "Low Nitrogen", "Pathogen")
#Now combining the new treatment columns into the old data frame
dat1 = cbind(treats$`Low Water`, treats$`Low Nitrogen`, treats$Pathogen, dat$Line, dat$Height_in)
head(dat1)
dat1 = as.data.frame(dat1)
dat1[,4] = sapply(dat1[,4],function(x){
  x = levels(dat$Line)[x]
})
head(dat1)
#give them column names
colnames(dat1) =c("Low Water", "Low Nitrogen", "Pathogen", "Line", "Height")
head(dat1)
#Write the csv
write.csv(dat1,file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/Manching BayesNet/ManchingScrubbed.csv",
          row.names = FALSE)

#####Remove any Mo017 or B73, also Mo066#####
length(unique(dat1$Line))
bad = which(dat1$Line == "Mo066" | 
              dat1$Line == "B73" |
              dat1$Line == "Mo17 parent" |
              dat1$Line == "IBMMo17" |
              dat1$Line == "IBMB73")
dat2 = dat1[-bad,]
length(unique(dat2$Line))
#then combine IMB lines with nonIBM lines of the same number
bad = which(substr(dat2$Line,1,3) == "IBM")
dat2 = dat2[-bad,]
length(unique(dat2$Line))

#getting rid of the space after Mo001
bad1 = which(dat2$Line == "Mo001 ")
dat2$Line[bad1] = "Mo001"
unique(dat2$Line)
#write the csv
write.csv(dat2,file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/Manching BayesNet/ManchingScrubbed.csv",
          row.names = FALSE)
#####Need to add in SNP info#####
snp = read.csv(file = "C:/Users/twili/Desktop/GIThub/Stapleton-Lab/vQTL IBM and Manching/Data/IBM94markerset08seq.csv")
relevant = data.frame(matrix(rep(0,length(dat2$Line)*dim(snp)[1]), ncol = dim(snp)[1]))
dat3 = sapply(dat2$Line, function(x){
  substr(x,2,2) = "O"
  column = which(colnames(snp) == x)
  vect = data.frame(as.character(snp[,column]))
  return(vect)
})
dat3 = as.data.frame(matrix(unlist(dat3), nrow = dim(dat2)[1], byrow = TRUE))
colnames(dat3) = colnames()
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