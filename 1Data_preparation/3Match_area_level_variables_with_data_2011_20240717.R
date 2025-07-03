######################################
##Script to match percentage greenspace dataset to area level variables
#The syntax was created by Zoe E Reed.
#The syntax was checked by Gareth J. Griffith
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################

library(rgeos)
library(sp)
library(rgdal)
library(raster)
library(data.table)
library(dplyr)

######################################
#Load data
######################################

#Read in MSOA shapefile and select relevant columns
m1<-readOGR("/path/Middle_Layer_Super_Output_Areas_(December_2011)_Boundaries/Middle_Layer_Super_Output_Areas_(December_2011)_Boundaries.shp", layer="Middle_Layer_Super_Output_Areas_(December_2011)_Boundaries")
m1_new<-m1[,-(c(1,3:6))]

#Read in IZ shapefile and select relevant columns
m2<-readOGR("/path/SG_IntermediateZoneBdry_2011/SG_IntermediateZone_Bdry_2011.shp", layer="SG_IntermediateZone_Bdry_2011")
m2_new<-m2[,-(c(2:9))]

#Change IZ column name to be same as that for MSOA
names(m2_new)<-"msoa11cd"

#Combine datasets
msoa<-raster::bind(m1_new, m2_new)

#Load UK biobank data
load("/path/green_bump_biobank_final_20240904.RData")

######################################
#Prepare data for matching to MSOA and LAD
######################################

data <- green_data

#Remove participants without east or north data
data<-data[!is.na(data$east),]
data<-data[!is.na(data$north),]

#Extract spatial coordinates of participants from UK biobank and set projection 
pts<-data.frame(data[,c("east", "north")])
pts2<-SpatialPoints(pts)
proj4string(pts2)<-CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=m +no_defs")

######################################
#Assign MSOA to participants
######################################

#Use overlap to see which regions people fall within
overlap<-sp::over(pts2, msoa)
sum(is.na(overlap$msoa11cd))

#Create a matrix the length of overlap for for loop below 
tmp<-matrix(nrow=dim(overlap)[1], ncol=1)

#For those with NA assign the neasrest MSOA
for(i in 1:length(tmp)){
	if(is.na(overlap$msoa11cd[i])){
		tmp[i]<-which.min(spDists(pts2[i], msoa))
	}
	else{
		tmp[i]<-NA
	}
}

for(i in 1:dim(overlap)[1]){
		if(is.na(overlap$msoa11cd[i])){
			overlap$msoa11cd[i]<-msoa$msoa11cd[as.numeric(tmp[i])]
			overlap$msoa11cd[i]<-msoa$msoa11cd[as.numeric(tmp[i])]
		}	
}

#Add column to main dataset for MSOA
data$census_msoa<-overlap$msoa11cd

######################################
#Save data for use in next script
######################################

save(data, file="/path/greenspace_biobank_final_census_20240904.RData")

