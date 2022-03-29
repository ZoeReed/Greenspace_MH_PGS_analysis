######################################
##Script to match NDVI data to main dataset and creates a dataset for NDVI
#The syntax was created by Zoe E Reed.
#The syntax was checked by ....
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

#Read in NDVI dataset
green<-read.csv("/filepath/UKBUMP_Greenness.csv", header=T, sep=",")
colnames(green)<-c("ID", colnames(green)[2:length(colnames(green))])

#Read in linker file to link UKBUMP to main phenotypic data
linker1<-read.csv("/filepath/Bridge_21829_ukbump_run11466.csv", header=T, sep=",")

#Read in main UKBB dataset from script 1
load("/filepath/phen_data_matched_distance_20200805.RData")

#Read in linker file to link genetic IDs data to main phenotypic data
linker2<-read.csv("/filepath/linker_ukb21829.csv", header=T, sep=",")
colnames(linker2)<-c("IID", "phen_ID")

#Merge all datasets to create one with NDVI, genetic IDs and phenotypic data
green_tmp<-merge(green, linker1, by.x="ID", by.y="eid_ukbump", all.x=T)
colnames(green_tmp)<-c(colnames(green), "phen_ID")
green_tmp2<-right_join(phen, green_tmp, by="phen_ID")
green_data<-left_join(green_tmp2, linker2, by="phen_ID")

######################################
#Tidy up dataset and NDVI variables
######################################

#Remove variables not needed
green_data$ID<-NULL
green_data$phen_ID<-NULL
green_data<-green_data[, c(18, 1:17)]
green_data$FID<-green_data$IID
green_data<-green_data[, c(19, 1:18)]

#Exclude participants with no IID (genetic data)
green_data<-green_data[!is.na(green_data$IID),]

#Where values are very small or large make these NA
green_data<-green_data[!(green_data$NDVI_500m_mean<(-1.800000e+100)) | is.na(green_data$NDVI_500m_mean),]
green_data<-green_data[!(green_data$NDVI_500m_min<(-1.800000e+100))| is.na(green_data$NDVI_500m_min),]
green_data<-green_data[!(green_data$NDVI_500m_max<(-1.800000e+100))| is.na(green_data$NDVI_500m_max),]
green_data<-green_data[!(green_data$NDVI_500m_STD<(-1.800000e+100))| is.na(green_data$NDVI_500m_STD),]
green_data<-green_data[!(green_data$NDVI_1000m_mean<(-1.800000e+100))| is.na(green_data$NDVI_1000m_mean),]
green_data<-green_data[!(green_data$NDVI_1000m_min<(-1.800000e+100))| is.na(green_data$NDVI_1000m_min),]
green_data<-green_data[!(green_data$NDVI_1000m_max<(-1.800000e+100))| is.na(green_data$NDVI_1000m_max),]
green_data<-green_data[!(green_data$NDVI_1000m_STD<(-1.800000e+100))| is.na(green_data$NDVI_1000m_STD),]

#Remove min 1000m as all -1 (not required for analyses anyway)
green_data$NDVI_1000m_min<-NULL

######################################
#Save data for use in next script
######################################

save(green_data, file="/filepath/green_bump_biobank_final_20200805.RData")

