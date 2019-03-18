### Graphs for height corn data ###
library("ggplot2")
library("car")

#setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing")
dat = read.csv(file = "dat2.csv", header = TRUE)
dat = dat[,c(3,4,7)]


### plot based on height and Mo###
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
dat.num = substrRight(as.integer(dat$Genotype), 3)
dat = cbind(dat.num, dat)

fullScat = scatterplot(dat$dat.num, dat$Height..in.., main="Scatterplot of Height based on Genotype", 
     xlab='Genotype', ylab="Height in inches", pch=19)

inbred = dat[1:896,]
hybrid = dat[-(1:896),]
par(mfrow = c(2,1))
InScat = plot(inbred$dat.num, inbred$Height..in.., main="Scatterplot of Height based on Genotype", 
              xlab='Genotype', ylab="Height in inches", pch=19)
OutScat = plot(hybrid$dat.num, hybrid$Height..in.., main="Scatterplot of Height based on Genotype", 
               xlab='Genotype', ylab="Height in inches", pch=19)

dev.off()

#######################
#####alternatively#####
#######################
dotchart(dat$Height..in.., main = "Dot Plot of Height")
par(mfrow = c(1,2))
dotchart(inbred$Height..in.., main = "Dot Plot of Inbred Height")
dotchart(hybrid$Height..in.., main = "Dot Plot of hybrid Height")
dev.off()

### histogram of all
hist(dat$Height..in.., main = "Histogram of Heights", col = "blue", xlab = "Height")
### plot based on inbred/hybrid 

hist(inbred$Height..in.., col=rgb(1,0,0,0.5),
     main="Overlapping Histogram of Inbred and hybrid Heights", xlab="Height")
hist(hybrid$Height..in.., col=rgb(0,0,1,0.5), add=T)
box()
legend("topleft", 
       c("Inbred", "hybrid"), 
       fill =c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)), 
       bty = "n")



#comparative box plot b/t inbred and hybrid

boxplot(split(dat$Height..in.., dat$BreedType), xlab = "BreedType", ylab="Plant Height",
        col = c("green", "blue"))


### summary stats of height for inbred and hybrid ###

summary(inbred$Height..in..)
summary(hybrid$Height..in..)

### Creating sample snpHeights ###
dat = read.csv(file = "snpHeight.csv", header = TRUE)


sampSnpHeight = dat[c(1:120, 896:986), 1:100]

write.csv(sampSnpHeight, file = "sampSnpHeight.csv")
