######################################
##Script to extract variables needed from main UKBB dataset for our project
#The syntax was created by Zoe E Reed.
#The syntax was checked by Gareth J. Griffith
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################
library(data.table)
library(plyr)
library(dplyr)
library(sp)

######################################
#Read in UKBB dataset
######################################
file<-fread("/filepath/data.11393.csv", data.table=F, na.strings="")

######################################
#Subset for variables of interest
######################################
vars<-c("eid", "54-0.0", "21003-0.0", "22001-0.0","22702-0.0", "22704-0.0")
phen<-file[,vars]

#Rename variables
colnames(phen)<-c("phen_ID","centre", "age", "sex", "east", "north")

#Change data type
phen$age<-as.numeric(phen$age)
phen$phen_ID<-as.integer(phen$phen_ID)

######################################
#Exclusions
######################################
#Remove those without ID
phen<-phen[!is.na(phen$phen_ID),]

#Remove those who have withdrawn consent
withdrawn<-read.csv("/filepath/w21829_20200204.csv", header=F)
phen<-phen[!(phen$phen_ID %in% withdrawn$V1), ]

######################################
#Match to assessment centre data
######################################

#Read in centre coordinates
centre<-read.csv("/filepath/AssessmentCentre_rounded_grid_coordinate.csv", header=T, sep=",")

#Match to centre for coordinates
test<-merge(phen, centre, by.x="centre", by.y="Assessment.centre", all.x=T)
phen<-test

#Obtain distance to centre (no used for these analyses)
phen$distance_east<-sqrt((as.numeric(phen$east) - as.numeric(phen$Rounded.grid.co.ordinate..Eastings))^2)
phen$distance_north<-sqrt((as.numeric(phen$north) - as.numeric(phen$Rounded.grid.co.ordinate..Northings))^2)
phen$distance<-sqrt((phen$distance_north)^2+(phen$distance_east)^2)
phen$distance_east<-NULL
phen$distance_north<-NULL
phen$distance_sqd<-phen$distance^2
phen<-phen[!is.na(phen$distance_sqd),]

######################################
#Save data for use in next script
######################################

save(phen, file="/filepath/phen_data_matched_distance_20200805.RData")

