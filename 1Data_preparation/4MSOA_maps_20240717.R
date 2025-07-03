######################################
##Script to create maps for MSOAs (or equivalent in Scotland)
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
library(RColorBrewer)

######################################
#Load data
######################################

load("/path/greenspace_biobank_final_census_20240904.RData")

msoa<-readOGR("/path/Census_MSOA.shp", layer="Census_MSOA")

######################################
#MSOA map
######################################

#Obtain list of populated MSOAs in UKBB and make data frame
msoa_populated<-unique(data$census_msoa)
test<-data.frame(table(data$census_msoa))

#Merge this data frame with the MSOA dataset
msoa_merged<-merge(msoa, test, by.x="msoa11cd", by.y="Var1")

#For those MSOAs not populated assign 0
msoa_merged$Freq[is.na(msoa_merged$Freq)]=0

#Create colours for mapping
test<-colorRamp(c("lightblue1", "mediumpurple"))
ramp.list<-rgb(test(seq(0,1, length=8)), max=255)
pal<-c("#FFFFFF", ramp.list)

tmp1<-cut(msoa_merged$Freq, breaks=c(0, 0.1, 100, 200, 300, 400, 500, 600, 700, max(msoa_merged$Freq)), include.lowest=T)

tmp2<-cut(msoa_merged$Freq, breaks=c(0, 0.1, 100, 200, 300, 400, 500, 600, 700, max(msoa_merged$Freq)), labels=pal, include.lowest=T)

msoa_merged$colour_index<-tmp1
msoa_merged$colours<-tmp2

#Create plot of UKBB population by MSOA and save as png file
png("/path/census_MSOA_20250114.png", width=1500, height=2000)

plot(msoa_merged, main="UK Biobank population by MSOA", col=pal[msoa_merged$colours], lwd=0.5, cex.main=3)

legend(470000, 1000000, legend=c("0", "1 to 100", "101 to 200", "201 to 300", "301 to 400", "401 to 500", "501 to 600", "601 to 700", "700+"), col=pal, pch=15, cex=2.3, pt.cex=2.3, title="N partcipants per MSOA")

dev.off()

