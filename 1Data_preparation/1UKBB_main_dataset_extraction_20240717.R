######################################
##Script to extract variables needed from main UKBB dataset for our project
#The syntax was created by Zoe E Reed.
#The syntax was checked by Gareth J. Griffith
#R Version 4.3.3
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

file<-fread("/path/data.51913.csv", data.table=F, na.strings="")

######################################
#Subset for variables of interest
######################################

vars<-c("eid", "54-0.0", "21003-0.0", "22001-0.0","22702-0.0", "22704-0.0", "24503-0.0")
phen<-file[,vars]

#Rename variables
colnames(phen)<-c("phen_ID","centre", "age", "sex", "east", "north", "greenspace_percent_300m")

#Change data type
phen$age<-as.numeric(phen$age)
phen$phen_ID<-as.integer(phen$phen_ID)

#Load linker file for genetic IDs
linker<-read.csv("/path/linker.81499.csv", header=T, sep=",")
head(linker)
colnames(linker)<-c("IID", "phen_ID")

#Match genetic data IDs to phenotypic data IDs using linker data
data<-merge(phen, linker, by="phen_ID", all.x=T)

######################################
#Exclusions
######################################

#Remove those without ID
data<-data[!is.na(phen$phen_ID),]

#Remove those withdrawn by date of creation
withdrawn<-read.csv("/path/w81499_2023-04-25.csv", header=F)
data<-data[!(data$phen_ID %in% withdrawn$V1), ]

######################################
#Match to assessment centre data
######################################

#Read in centre coordinates
centre<-read.csv("/path/AssessmentCentre_rounded_grid_coordinate.csv", header=T, sep=",")

#Match to centre for coordinates
tmp<-merge(data, centre, by.x="centre", by.y="Assessment.centre", all.x=T)
data<-tmp

#Obtain distance to centre (no used for these analyses)
data$distance_east<-sqrt((as.numeric(data$east) - as.numeric(data$Rounded.grid.co.ordinate..Eastings))^2)
data$distance_north<-sqrt((as.numeric(data$north) - as.numeric(data$Rounded.grid.co.ordinate..Northings))^2)
data$distance<-sqrt((data$distance_north)^2+(data$distance_east)^2)
data$distance_east<-NULL
data$distance_north<-NULL
data$distance_sqd<-data$distance^2
data<-data[!is.na(data$distance_sqd),]

######################################
#Format data
######################################

data<-data[!is.na(data$IID),]

######################################
#Save data for use in next script
######################################

save(data, file="/path/greenspace_biobank_final_20240904.RData")

