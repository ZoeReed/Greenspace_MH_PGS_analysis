######################################
##Script to estimate R squared for PCs and plot these
#The syntax was created by Zoe E Reed.
#The syntax was checked by Gareth J. Griffith
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################

library(plyr)
library(ggplot2)

######################################
#Load data
######################################

#NDVI
load("/path/greenspace_biobank_final_census_20240904.RData")
ndvi<-data
rm(data)

#Principal components (PCs)
pcs<-read.table(file="/path/100PCs_covariates.txt", header=T)

#Match data
common.samples<-intersect(as.character(ndvi$IID), as.character(pcs$IID))
ndvi<-ndvi[match(common.samples, ndvi$IID),]
pcs<-pcs[match(common.samples, pcs$IID),]

######################################
#Create models for NDVI and PCs
######################################

#Create empty list for results
ndvi_results<-list()

#Run first model for PC1
model<-lm(ndvi$NDVI_500m_mean~pcs$PC1)

#Make data frame from this data to be filled in loop below for remaining PCs
ndvi_tmp<-as.data.frame(1)
ndvi_tmp$ndvi_rsq<-summary(model)$r.squared
colnames(ndvi_tmp)<-c("PC", "ndvi_rsq")

#For loop for models with PC2 to PC100
for(i in 4:102){
	print(i-2)
	model<-update(model, as.formula(paste0(".~. + pcs[,",i,"]")))
	new<-data.frame(matrix(, nrow=1, ncol=0))
	new$PC<-i-2
	new$ndvi_rsq<-summary(model)$r.squared
	ndvi_results<-rbind(ndvi_results, new)
}

#Combine into one results dataset
ndvi_results_final<-rbind(ndvi_tmp, ndvi_results)

######################################
#Save data
######################################

write.table(ndvi_results_final, file ="/path/ndvi_pc_R2.csv", row.names = FALSE, append = FALSE, col.names = TRUE, sep = "\t ", quote=F)

######################################
#Load data for Percentage greenspace
######################################

#Percentage greenspace
load("/path/greenspace_biobank_final_census_20240904.RData")
greenspace<-data
rm(data)

#PCs
pcs<-read.table(file="/path/100PCs_covariates.txt", header=T)

#Match data
common.samples<-intersect(as.character(greenspace$IID), as.character(pcs$IID))
greenspace<-greenspace[match(common.samples, greenspace$IID),]
pcs<-pcs[match(common.samples, pcs$IID),]

######################################
#Create models for NDVI and PCs
######################################

#Create empty list for results
greenspace_results<-list()

#Run first model for PC1
model<-lm(greenspace$greenspace_percent_300m~pcs$PC1)

#Make data frame from this data to be filled in loop below for remaining PCs
greenspace_tmp<-as.data.frame(1)
greenspace_tmp$greenspace_rsq<-summary(model)$r.squared
colnames(greenspace_tmp)<-c("PC", "greenspace_rsq")

rm(i)

#For loop for models with PC2 to PC100
for(i in 4:102){
	model<-update(model, as.formula(paste0(".~. + pcs[,",i,"]")))
	new<-data.frame(matrix(, nrow=1, ncol=0))
	new$PC<-i-2
	new$greenspace_rsq<-summary(model)$r.squared
	greenspace_results<-rbind(greenspace_results, new)
}

#Combine into one results dataset
greenspace_results_final<-rbind(greenspace_tmp, greenspace_results)

######################################
#Save data
######################################

write.table(greenspace_results_final, file ="/path/greenspace_pc_R2.csv", row.names = FALSE, append = FALSE, col.names = TRUE, sep = "\t ", quote=F)

######################################
#Create plot for both NDVI and percentage greenspace
######################################

#Merge datasets
all<-join_all(list(ndvi_results_final, greenspace_results_final), by='PC', type='full')

#Create plot and save as a pdf
points<-seq(0,100, 5)	
	
pdf("/path/rsquared_for_pcs_greenness.pdf",width=20, height=10)

ggplot(data=all, aes(x=PC)) +
	theme(text = element_text(size=25)) +
	geom_line(aes(y=ndvi_rsq, colour="NDVI"), size=1.1) +
	geom_line(aes(y=greenspace_rsq, colour="Percentage greenspace"), size=1.1) +
	ylab(label="R squared") +
	xlab("Number of Principal Components") +
	scale_x_continuous(breaks=points, expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0)) +
	scale_colour_manual("", values=c("NDVI"="azure3", "Percentage greenspace"="plum3"))
	
dev.off()	
	