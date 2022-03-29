######################################
##Script to match percentage greenspace data to main dataset and creates a dataset for percentage greenspace
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

#Read in percentage greenspace dataset
green<-read.csv("/filepath/data.37201.csv", header=T, sep=",", na.strings="")

#Rename columns
colnames(green)<-c("phen_ID","greenspace_percent_1000m", "greenspace_percent_1000m_2nd", "Domestic_garden_percent_1000m", "Domestic_garden_percent_1000m_2nd", "Water_percent_1000m", "Water_percent_1000m_2nd", "greenspace_percent_300m", "greenspace_percent_300m_2nd", "Domestic_garden_percent_300m", "Domestic_garden_percent_300m_2nd", "Water_percent_300m", "Water_percent_300m_2nd", "natural_env_percent_1000m", "natural_env_percent_1000m_2nd", "natural_env_percent_300m", "natural_env_percent_300m_2nd", "distance_to_coast", "distance_to_coast_2nd", "IMD_england", "IMD_wales", "IMD_scotland")


#Read in main UKBB dataset from script 1
load("/filepath/phen_data_matched_distance_20200805.RData")

#Read in linker file to link genetic IDs data to main phenotypic data
linker2<-read.csv("/filepath/linker_ukb21829.csv", header=T, sep=",")
colnames(linker2)<-c("IID", "phen_ID")

#Merge all datasets to create one with percentage greenspace, genetic IDs and phenotypic data
green_tmp<-right_join(phen, green, by="phen_ID")
green_data<-left_join(green_tmp, linker2, by="phen_ID")

#make sure all withdrawn
withdrawn<-read.csv("/filepath/w21829_20200204.csv", header=F)
green_data<-green_data[!(green_data$phen_ID %in% withdrawn$V1), ]

######################################
#Tidy up dataset and percentage greenspace variables
######################################

#Remove variables not needed
green_data$phen_ID<-NULL
green_data<-green_data[, c(32, 1:31)]
green_data$FID<-green_data$IID
green_data<-green_data[, c(33, 1:32)]

#Exclude participants with no IID (genetic data)
green_data<-green_data[!is.na(green_data$IID),]

######################################
#Save data for use in next script
######################################

save(green_data, file="/filepath/greenspace_biobank_final_20200805.RData")

