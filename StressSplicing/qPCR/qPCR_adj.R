### QPCR PLATE FRAMING & ADJUSTMENT MODEL ###
library(tidyr)
library(pracma)
library(stringr)
library(tidyverse)
library(dplyr)
# Mac Directory
#setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_11")
#setwd("~/Stapleton_Lab/Stapleton_Lab/Stress_Splicing/2018_(MONTH)")
# PC Directory
#setwd(~/Desktop/GIThub/StapletonLab/StressSplicing/qPCR/)

### READ IN DERIVATIVE DATA###
# In the case of having two separate CSV files of calculated derivatives,
# use this code to combine, prior to the following transpositions:
deriv.1<-read.csv(file = "2018_11_1_plate_qPCR_output.csv", header=FALSE)
deriv.2<-read.csv(file = "2018_11_2_plate_qPCR_output.csv", header=FALSE)
deriv=cbind(deriv.1, deriv.2)

# In the case of having one CSV containing calculated derivatives, use this code:
#deriv=read.csv(file = "(YEAR_MONTH_PLATE_qPCR_output.csv", header=FALSE)

########################################################## 
################### Initial Data Framing #################
########################################################## 

# Remove extra labels row and column 
deriv = deriv[-1,-1]
# Transpose derivatives to be in equivalent format as raw plate data
deriv = as.data.frame(t(deriv), header=TRUE)
# Rename columns
colnames(deriv)=c("reaction_type", "sampleID", "starting_quantity", "cpD1", "cpD2")
# Remove extra labels row
deriv=deriv[-1,]

############################
## Removing NTC and Minus ##
############################

# Indicate if sample is NTC (negative control)
deriv['sampleID_NTC'] = grepl('NTC', deriv$sampleID)
# Remove NTC samples, indicator (T/F) column, and cpD2 values
ntc = which(deriv$sampleID_NTC)
deriv = deriv[-ntc,]
deriv = deriv[,-c(5,6)]
# Indicate if sample is 'Plus' or 'Minus'
deriv['sampleID_Minus'] = grepl('minus', deriv$sampleID)
# Remove 'Minus' values (include only gblock+ values), and indicator (T/F) column
minus = which(deriv$sampleID_Minus)
deriv = deriv[-minus,]
deriv = deriv[,-c(5)]

### COMPLETED INITIAL DATA FRAMING ###

########################################################## 
################# Calibrated Data Framing ################
########################################################## 

# Create/Write data frame for Calibrated values
calib_df = deriv %>% filter(str_detect(sampleID, "g"))
# Sort by starting quantity
calib_df = calib_df[order(calib_df$starting_quantity),]
calib_df$starting_quantity = as.numeric(as.character(calib_df$starting_quantity))
calib_df$cpD1 = as.numeric(as.character(calib_df$cpD1))
calib_data = calib_df


# Create empty vectors for for-loop to input cpD1 values
test1 = c()
allP = c()
startq = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(calib_data$starting_quantity)){
  sq <- calib_data$starting_quantity[i]
  if(i %% 6 == 1){
    startq = c(startq,sq,sq,sq)
  }
  val <- toString(calib_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1 = c(test1, calib_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP = c(allP, calib_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by starting quantity
calib_data = cbind(startq, test1, allP)
# Format starting quantity values as decimals, not scientific notation
calib_data[,1]=format(calib_data[,1], scientific=FALSE)

### COMPLETED CALIBRATED DATA FRAME ###

########################################################## 
############ ADJUSTMENT MODEL Calibrated Data ############
########################################################## 

#data=read.csv(file="2018_11_Calibrated_Data_Frame.csv", header=TRUE) 
# Create empty vectors for for-loop input
data = as.data.frame(calib_data)
data$test1 = as.numeric(as.character(data$test1))
data$allP = as.numeric(as.character(data$allP))
adj_val = c()
allP = c()
startq = c()
ratio =data$allP/data$test1
# Itterating through each set of (3) observations performing U-Stats on each set of inputs
for (i in 1:(nrow(data)/3)){
  t_x <- c(data$allP[3*i - 2], data$allP[3*i - 1], data$allP[3*i])
  t_y <- c(data$test1[3*i - 2], data$test1[3*i - 1], data$test1[3*i])
  adj <- mean(outer(t_x, t_y, "-"))
  adj_val <- c(adj_val, adj, adj, adj)
}
adjusted_test1 <- test1 + adj_val 
# Creating the adjustment model lm(y-axis~x-axis)
# Changed adj_val^2 to adj_val to try to make the model better --> VERIFY THIS WITH DR. WANG
adj_model <- lm(adj_val~ratio) #Adjusted/avg slopes model --> to get JC VQTL vals 
summary(adj_model) 
par(mfrow = c(2,2))
plot(adj_model)
dev.off()
scatter.smooth(ratio, adj_val)
abline(ratio, adj_val)


compare_data = cbind(adj_val, ratio)
write.csv(file = "compare_dataJulia.csv", compare_data)
### COMPLETED ADJUSTMENT MODEL - CALIBRATED DATA ###

### EXPERIMENTAL DATA FRAME ###
# Create/Write data frame for Calibrated values
exp_df = deriv %>% filter(str_detect(sampleID, "g")==FALSE)
# Sort by starting quantity
exp_df = exp_df[order(exp_df$starting_quantity),]
# Remove first and last rows (unnecessary labeling)
exp_df = exp_df[-1,]
exp_df = exp_df[-nrow(exp_df),]
write.csv(exp_df, file="2018_11_Experimental_Data_Frame.csv")
# Create new transposed data set
exp_data = read.csv(file="2018_11_Experimental_Data_Frame.csv", header=TRUE) 
# Remove extra labeling column
exp_data = exp_data[,-1]
# Order data by sampleID
exp_data = exp_data[order(exp_data$sampleID),]
# expdata2 will be used to compare to the for-loop output, exp_data, to verify processed correctly
expdata2=exp_data

#REMOVE SAMPLES WITH MISSING ALLP/TEST1 VALUES; 
#CONFIRM WITH DR.S IF SAMPLES WITH UNUSUAL CP VALUES SHOULD BE OMITTED
exp_data=exp_data ##Finish this code to remove certain sampleID's -- PRIOR to for-loop

# Create empty vectors for for-loop to input cpD1 values
test1.exp = c()
allP.exp = c()
sampleID.exp = c()
# For loop -- iterating thru starting quantity and reaction type to return cpD1 values 
for(i in 1:length(exp_data$sampleID)){
  id.exp = toString(exp_data$sampleID[i])
  if(i %% 2 == 1){
    sampleID.exp = c(sampleID.exp, id.exp)
  }
  val = toString(exp_data$reaction_type[i])
  if(strcmp(val, "test1")){
    test1.exp = c(test1.exp, exp_data$cpD1[i])
  }
  if(strcmp(val, "all_products")){
    allP.exp = c(allP.exp, exp_data$cpD1[i])
  }
}
# Bind test1 and allProd cpD1 values by sample ID
exp_data = cbind(sampleID.exp, test1.exp, allP.exp)
# Write CSV file
write.csv(exp_data, file="2018_11_Experimental_Data_Frame.csv")
### COMPLETED EXPERIMENTAL DATA FRAME ###

### ADJUSTMENT MODEL - EXPERIMENTAL DATA ### 
# Using the adjustment model on the expiremental data
new = data.frame(ratio = exp_data$all_productsPrimers_Cp1/exp_data$test1_Cp1)
predict(adj_model, new , interval = "confidence")
#---> Fill in any remaining parts of experimental adjusmtent model

### COMPLETED ADJUSTMENT MODEL - EXPERIMENTAL DATA ###