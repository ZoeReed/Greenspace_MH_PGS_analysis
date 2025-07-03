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

library(stats)
library(DescTools)

######################################
#Read in data
######################################

#Phenotypic
load("/path/MH_biobank_final_20240904.RData")

#GRS
well_grs<-read.csv("/path/grs_well_05.csv")

######################################
#Format data
######################################

#Merge datasets into a single data frame
all<-merge(phen, well_grs, by.x="IID", by.y="id", all.x=T) 

dim(all)
head(all)
str(all)

phen <- all

phen$sex<-as.factor(phen$sex)

phen$wellbeing<-as.numeric(phen$wellbeing)

######################################
#Exclusions
######################################

withdrawn<-read.csv("/path/w81499_2023-04-25.csv", header=F)
phen<-phen[!(phen$eid %in% withdrawn$V1), ]

dim(phen)

rec_exc<-read.table("/path/data.combined_recommended.qctools.txt", header=F)
phen<-phen[!(phen$IID %in% rec_exc$V1), ]
dim(phen)

high_rel<-read.csv("/path/data.highly_relateds.qctools.txt", header=F)
phen<-phen[!(phen$IID %in% high_rel$V1), ]
dim(phen)

min_rel<-read.csv("/path/data.minimal_relateds.qctools.txt", header=F)
phen<-phen[!(phen$IID %in% min_rel$V1), ]
dim(phen)

non_white_brit<-read.csv("/path/data.non_white_british.qctools.txt", header=F)
phen<-phen[!(phen$IID %in% non_white_brit$V1), ]
dim(phen)

#Restrict to only those with genetic data
phen<-phen[!is.na(phen$grs),]
dim(phen)

######################################
#Standardise
######################################

#Standardise NDVI
phen$NDVI_500m_mean<-scale(phen$NDVI_500m_mean, center=T, scale=T)

#Standardise percentage greenspace
phen$greenspace_percent_300m<-scale(phen$greenspace_percent_300m, center=T, scale=T)

#Standardise PGS
phen$grs<-scale(phen$grs, center=T, scale=T)

######################################
#R-squared
######################################

#R-squared for wellbeing
mod<-lm(phen$wellbeing ~ phen$grs)
summary(mod)$r.squared

######################################
#Phenotypic analyses
######################################

#Association between depression diagnosis and NDVI, adjusted for sex and age
mod<-lm(phen$NDVI_500m_mean ~ phen$wellbeing + phen$sex + phen$age)
summary(mod)
nobs(mod)
confint.default(mod)

#Association between depression diagnosis and percentage greenspace, adjusted for sex and age
mod<-lm(phen$greenspace_percent_300m ~ phen$wellbeing + phen$sex + phen$age)
summary(mod)
nobs(mod)
confint.default(mod)

