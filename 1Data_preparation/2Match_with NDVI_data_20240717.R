######################################
##Script to match NDVI data to main dataset and creates a dataset for NDVI
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

######################################
#Read in datasets
######################################

#Read in NDVI dataset, there are several files and a seperate header file, use NDVI 500m mean
birm<-read.csv("/path/Files for RETMAN/Folder_5_Greenness/Folder_5_Greenness/UKB_Birmingham_Nottingham_NDVI.csv", header=T)

bris<-read.csv("/path/Built Environment - Bristol/Folder_5_Greenness/UKB_Bristol_NDVI.csv", header=T)

lon<-read.csv("/path/Files for Reman GLA/Folder_5_Greenness/UKB_London_NDVI.csv", header=T)
lon <- lon[,1:5]

LS<-read.csv("/path/Files for Retman - LS/Folder_5_Greenness/UKB_Leeds_Sheffield_NDVI.csv", header=T)

liv1<-read.csv("/path/Files for Retman - Liverpool Manchester Bury_New/Folder_5_Greenness/UKB_LIV_MAN_BURY_NDVI500m.csv", header=T)

Stoke<-read.csv("/path/Files for Retman - Stoke/Folder_5_Greenness/UKB_Stoke_NDVI500m.csv", header=T)

Wales<-read.csv("/path/Files for Retman - Wales/Files_5_Greenness/Wales_UKB_NDVI.csv", header=T)
Wales <- Wales[,1:5]

combined <- rbind(birm, bris, lon, LS, liv1, Stoke, Wales)

header<-read.csv("/path/Files for RETMAN/Folder_5_Greenness/Folder_5_Greenness/UKB_Birmingham_Nottingham_NDVI_Header.csv", header=T)

colnames(combined)<-c(header)

combined <- combined[,c("Encoded anonymised participant ID", "NDVI_500m_mean")]

colnames(combined)<-c("bumpID", "NDVI_500m_mean")

#Read in main UKBB dataset from script 1
load("/path/greenspace_biobank_final.RData")

#Read in linker file to link UKBUMP to main phenotypic data
linker1<-read.table("/path/ukb81499bridge38006.txt", header=F, sep=" ")
colnames(linker1)<-c("phen_ID", "bumpID")

#Merge all datasets to create one with NDVI, genetic IDs and phenotypic data
green_tmp<-merge(combined, linker1, by="bumpID", all.x=T)
green_data<-merge(data, green_tmp, by="phen_ID", all.x=T)

######################################
#Tidy up dataset and NDVI variables
######################################

#Remove variables not needed
green_data$ID<-NULL
green_data$phen_ID<-NULL
green_data<-green_data[, c(7, 1:6, 13)]
green_data$FID<-green_data$IID
green_data<-green_data[, c(9, 1:8)]

#Where values are very small or large make these NA
green_data<-green_data[!(green_data$NDVI_500m_mean<(-1.800000e+100)) | is.na(green_data$NDVI_500m_mean),]

######################################
#Save data for use in next script
######################################

save(green_data, file="/path/green_bump_biobank_final_20240904.RData")

