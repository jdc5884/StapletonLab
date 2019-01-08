### qPCR Analysis on ER Stress splicing ###

setwd("C:/Users/twili/Desktop/GIThub/StapletonLab/StressSplicing")
Plate = read.csv(file = "2018_6_1_plate.csv", header = TRUE)

Cycle.m = Plate[6:40,]
