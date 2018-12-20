rcsv<- read.csv(file = url("https://raw.githubusercontent.com/tbillman/Stapleton-Lab/master/vQTL%20Random%20and%20Family/data/tidied/Random2.csv"))
fcsv<- read.csv(file = url("https://raw.githubusercontent.com/tbillman/Stapleton-Lab/master/vQTL%20Random%20and%20Family/data/tidied/Family.csv"))
fcsv <- fcsv[-1:-2,]
PlantHeight <- c(rcsv$height.in.,fcsv$PlantHeight)
TasselBranchNum <- c(rcsv$NumTasselBranches, fcsv$tasselbranchNum)
TasselBranchDeg <- c(rcsv$TasselBranchAngle.deg., fcsv$tasselbranchAngle_degrees)
gdat <- rbind(rcsv[,-1:-3], fcsv[,-1:-3])
RFcovar <- c("NA","NA", rep(0,dim(rcsv)[1]-2), rep(1,dim(fcsv)[1]))
fdat <- cbind(PlantHeight, TasselBranchNum, TasselBranchDeg, RFcovar, gdat)
write.csv(fdat, file = "C://Users/Thomas/Documents/GitHub/Stapleton-Lab/vQTL\ Random\ and\ Family/combined/combined.csv",
          row.names = F)
