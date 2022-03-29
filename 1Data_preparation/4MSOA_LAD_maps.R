######################################
##Script to create maps for MSOA and LAD (or equivalent in Scotland)
#The syntax was created by Zoe E Reed.
#The syntax was checked by Gareth J. Griffith
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################
library(RColorBrewer)

######################################
#Load data
######################################
load("/filepath/greenspace_biobank_final_census_20200805.RData")

msoa<-readOGR("/filepath/Census_MSOA.shp", layer="Census_MSOA")

l1<-readOGR("/filepath/Local_Authority_Districts__December_2017__Boundaries_in_Great_Britain.shp", layer="Local_Authority_Districts__December_2017__Boundaries_in_Great_Britain")
lad<-l1[,-(c(1,3:8))]


######################################
#MSOA map
######################################

#Obtain list of populated MSOAs in UKBB and make data frame
msoa_populated<-unique(green_data$census_msoa)
test<-data.frame(table(green_data$census_msoa))

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
png("/filepath/census_MSOA.png", width=1500, height=2000)

plot(msoa_merged, main="UK Biobank population by MSOA", col=pal[msoa_merged$colours], lwd=0.5, cex.main=3)

legend(470000, 1000000, legend=c("0", "1 to 100", "101 to 200", "201 to 300", "301 to 400", "401 to 500", "501 to 600", "601 to 700", "700+"), col=pal, pch=15, cex=2.3, pt.cex=2.3, title="N partcipants per MSOA")

dev.off()


######################################
#LAD map
######################################

#Remove datasets for MSOA
rm(test)
rm(tmp1)
rm(tmp2)

#Obtain list of populated LADs in UKBB and make data frame
lad_populated<-unique(green_data$census_lad)
test<-data.frame(table(green_data$census_lad))

#Merge this data frame with the LAD dataset
lad_merged<-merge(lad, test, by.x="lad17cd", by.y="Var1")

#For those LADs not populated assign 0
lad_merged$Freq[is.na(lad_merged$Freq)]=0

#Create colours for mapping
tmp1<-cut(lad_merged$Freq, breaks=c(0, 0.1, 2500, 5000, 7500, 10000, 12500, 15000, 17500, max(lad_merged$Freq)), include.lowest=T)

tmp2<-cut(lad_merged$Freq, breaks=c(0, 0.1, 2500, 5000, 7500, 10000, 12500, 15000, 17500, max(lad_merged$Freq)), labels=pal, include.lowest=T)

lad_merged$colour_index<-tmp1
lad_merged$colours<-tmp2

#Create plot of UKBB population by LAD and save as png file
png("/filepath/census_LAD.png", width=1500, height=2000)

plot(lad_merged, main="UK Biobank population by LAD", col=pal[lad_merged$colours], lwd=0.5, cex.main=3)

legend(470000, 1000000, legend=c("0", "1 to 2,500", "2,501 to 5,000", "5,001 to 7,500", "7,501 to 10,000", "10,001 to 12,500", "12,501 to 15,000", "15,001 to 17,500", "17,500+"), col=pal, pch=15, cex=2.3, pt.cex=2.3, title="N partcipants per LAD")

dev.off()

