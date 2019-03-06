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
outbred = dat[-(1:896),]
par(mfrow = c(2,1))
InScat = plot(inbred$dat.num, inbred$Height..in.., main="Scatterplot of Height based on Genotype", 
              xlab='Genotype', ylab="Height in inches", pch=19)
OutScat = plot(outbred$dat.num, outbred$Height..in.., main="Scatterplot of Height based on Genotype", 
               xlab='Genotype', ylab="Height in inches", pch=19)

dev.off()

#######################
#####alternatively#####
#######################
dotchart(dat$Height..in.., main = "Dot Plot of Height")
par(mfrow = c(1,2))
dotchart(inbred$Height..in.., main = "Dot Plot of Inbred Height")
dotchart(outbred$Height..in.., main = "Dot Plot of Outbred Height")
dev.off()

### histogram of all
hist(dat$Height..in.., main = "Histogram of Heights", col = "blue", xlab = "Height")
### plot based on inbred/outbred 
inbred = dat[1:896,]
outbred = dat[-(1:896),]

hist(inbred$Height..in.., col=rgb(1,0,0,0.5),
     main="Overlapping Histogram of Inbred and Outbred Heights", xlab="Height")
hist(outbred$Height..in.., col=rgb(0,0,1,0.5), add=T)
box()
legend("topleft", 
       c("Inbred", "Outbred"), 
       fill =c(rgb(1,0,0,0.5),rgb(0,0,1,0.5)), 
       bty = "n")



#comparative box plot
values = merge(inbred, outbred, by = "")



