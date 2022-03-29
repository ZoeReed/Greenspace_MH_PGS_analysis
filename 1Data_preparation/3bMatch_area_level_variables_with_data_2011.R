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
library(rgdal) #add
library(raster)
library(data.table)
library(dplyr)

######################################
#Load data
######################################

#Read in MSOA shapefile and select relevant columns
m1<-readOGR("/filepath/Middle_Layer_Super_Output_Areas_(December_2011)_Boundaries.shp", layer="Middle_Layer_Super_Output_Areas_(December_2011)_Boundaries")
m1_new<-m1[,-(c(1,3:6))]

#Read in IZ shapefile and select relevant columns
m2<-readOGR("/filepath/SG_IntermediateZone_Bdry_2011.shp", layer="SG_IntermediateZone_Bdry_2011")
m2_new<-m2[,-(c(2:9))]

#Change IZ column name to be same as that for MSOA
names(m2_new)<-"msoa11cd"

#Combine datasets
msoa<-raster::bind(m1_new, m2_new)

#Read LAD file and select relevant columns
l1<-readOGR("/filepath/Local_Authority_Districts__December_2017__Boundaries_in_Great_Britain.shp", layer="Local_Authority_Districts__December_2017__Boundaries_in_Great_Britain")
lad<-l1[,-(c(1,3:8))]

#Load UK biobank data for percentage greenspace (from script 2a)
load("/filepath/greenspace_biobank_final_20200805.RData")

######################################
#Prepare data for matching to MSOA and LAD
######################################

#Remove participants without east or north data
green_data<-green_data[!is.na(green_data$east),]
green_data<-green_data[!is.na(green_data$north),]

#Extract spatial coordinates of participants from UK biobank and set projection 
pts<-data.frame(green_data[,6:7])
pts2<-SpatialPoints(pts)
proj4string(pts2)<-CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs")


######################################
#Assign MSOA to participants
######################################

#Use overlap to see which regions people fall within
overlap<-sp::over(pts2, msoa)
sum(is.na(overlap$msoa11cd)) #number of those who fall outside areas = 

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
green_data$census_msoa<-overlap$msoa11cd

######################################
#Assign LAD to participants
######################################

#Load in OA lookup table to match MSOA to LAD
lookup<-read.csv("/filepath/Output_Area_to_LSOA_to_MSOA_to_Local_Authority_District_(December_2017)_Lookup_with_Area_Classifications_in_Great_Britain.csv", header=T)

lookup2<-lookup[!duplicated(lookup$MSOA11CD),]

dim(green_data)
head(green_data)

#Merge datasets
green_data2<-merge(green_data, lookup2[,c("MSOA11CD", "LAD17CD")], by.x="census_msoa", by.y="MSOA11CD", all.x=T)

dim(green_data2)
head(green_data2)

#Rename column for LAD
colnames(green_data2)[colnames(green_data2) == 'LAD17CD']<-'census_lad'

green_data<-green_data2

######################################
#Save data for use in next script
######################################

save(green_data, file="/filepath/greenspace_biobank_final_census_20200805.RData")

