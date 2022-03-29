######################################
######################################
##Script to obtain R squared for wellbeing PGS and phenotypic analyses
#The syntax was created by Zoe E Reed.
#The syntax was checked by Tim T. Morris
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################
library(AER)
library(stats)

######################################
#Read in data
######################################

#Phenotypic
load("/filepath/phenotypic_one_sample_MR.RData")

#PGS
well_grs<-read.csv("/filepath/grs_well_05.csv")

######################################
#Format data
######################################

#Merge datasets into a single data frame
all<-merge(phen, well_grs, by.x="IID", by.y="id") 

all$sex<-as.factor(all$sex)

#treat wellbeing as continuous phenotype
all$wellbeing<-as.numeric(all$wellbeing)

######################################
#Load PCs
######################################

pcs_all<-read.table("/filepath/data.pca1-40.plink.txt")
colnames(pcs_all)<-c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", "PC33", "PC34", "PC35", "PC36", "PC37", "PC38", "PC39", "PC40")
vars<-c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25")
pcs<-pcs_all[,vars]

#Merge data
tmp<-all
all<-merge(tmp, pcs, by="IID", all.x=T)

######################################
#Exclusions
######################################

#Make sure all withdrawn are removed
withdrawn<-read.csv("/filepath/w21829_20200204.csv", header=F)
all<-all[!(all$eid %in% withdrawn$V1), ]

dim(all)

rec_exc<-read.table("/filepath/data.combined_recommended.qctools.txt", header=F)
all<-all[!(all$IID %in% rec_exc$V1), ]
dim(all)

high_rel<-read.csv("/filepath/data.highly_relateds.qctools.txt", header=F)
all<-all[!(all$IID %in% high_rel$V1), ]
dim(all)

min_rel<-read.csv("/filepath/data.minimal_relateds.qctools.txt", header=F)
all<-all[!(all$IID %in% min_rel$V1), ]
dim(all)

non_white_brit<-read.csv("/filepath/data.non_white_british.qctools.txt", header=F)
all<-all[!(all$IID %in% non_white_brit$V1), ]
dim(all)

#Restrict to only those with genetic data
all<-all[!is.na(all$grs),]
dim(all) #336,997

#Restrict to only those with genetic data
all<-all[!is.na(all$grs),]
dim(all) #336,997

######################################
#Standardise
######################################

#Standardise NDVI
all$NDVI_500m_mean<-scale(all$NDVI_500m_mean, center=T, scale=T)

#Standardise percentage greenspace
all$greenspace_percent_300m<-scale(all$greenspace_percent_300m, center=T, scale=T)

#Standardise PGS
all$grs<-scale(all$grs, center=T, scale=T)

######################################
#R-squared
######################################

#R-squared for wellbeing
mod<-lm(all$wellbeing ~ all$grs)
summary(mod)$r.squared

######################################
#Phenotypic analyses
######################################

#Association between depression diagnosis and NDVI, adjusted for sex and age
mod<-lm(all$NDVI_500m_mean ~ all$wellbeing + all2$sex + all2$age)
summary(mod)
nobs(mod)

#Association between depression diagnosis and percentage greenspace, adjusted for sex and age
mod<-lm(all$greenspace_percent_300m ~ all$wellbeing + all2$sex + all2$age)
summary(mod)
nobs(mod)
