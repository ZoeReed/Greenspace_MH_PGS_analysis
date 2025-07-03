######################################
##Script to obtain R squared for depression PGS and phenotypic analyses
#The syntax was created by Zoe E Reed.
#The syntax was checked by Tim T. Morris
#R Version 3.6.2
#The code below shows phen manipulations and recodes with annotations

######################################
#Load libraries (use instphen.packages() to instphen)
######################################

library(stats)
library(DescTools)

######################################
#Read in data
######################################

#Phenotypic
load("/path/MH_biobank_final_20240904.RData")

#GRS
dep_grs<-read.csv("/path/grs_dep_05.csv")

######################################
#Format data
######################################

all<-merge(phen, dep_grs, by.x="IID", by.y="id", all.x=T) 

dim(all)
head(all)
str(all)

phen <- all

phen$sex<-as.factor(phen$sex)
phen$MHQ_dep<-as.factor(phen$MHQ_dep)
phen$depression_diagnosis<-as.factor(phen$depression_diagnosis)

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
#Descriptives for all data
######################################

sum(!is.na(phen$depression_diagnosis))
table(phen$depression_diagnosis)

sum(!is.na(phen$MHQ_dep))
table(phen$MHQ_dep)

sum(!is.na(phen$wellbeing))
mean(as.numeric(phen$wellbeing), na.rm=T)
sd(as.numeric(phen$wellbeing), na.rm=T)

sum(!is.na(phen$schizophrenia_diagnosis))
table(phen$schizophrenia_diagnosis)

sum(!is.na(phen$age))
mean(as.numeric(phen$age), na.rm=T)
sd(as.numeric(phen$age), na.rm=T)

sum(!is.na(phen$sex))
table(phen$sex)

sum(!is.na(phen$NDVI_500m_mean))
mean(as.numeric(phen$NDVI_500m_mean), na.rm=T)
sd(as.numeric(phen$NDVI_500m_mean), na.rm=T)
range(phen$NDVI_500m_mean, na.rm=T)

sum(!is.na(phen$greenspace_percent_300m))
mean(as.numeric(phen$greenspace_percent_300m), na.rm=T)
sd(as.numeric(phen$greenspace_percent_300m), na.rm=T)

cor.test(phen$NDVI_500m_mean, phen$greenspace_percent_300m)

######################################
#Plot where people with NDVI live
######################################

library(rgeos)
library(sp)
library(rgdal)
library(raster)
library(RColorBrewer)

msoa<-readOGR("/path/Census_MSOA.shp", layer="Census_MSOA")

l1<-readOGR("/path/Local_Authority_Districts__December_2017__Boundaries_in_Great_Britain.shp", layer="Local_Authority_Districts__December_2017__Boundaries_in_Great_Britain")
lad<-l1[,-(c(1,3:8))]

tmp <- phen[!is.na(phen$NDVI_500m_mean),]

msoa_populated<-unique(tmp$census_msoa)
test<-data.frame(table(tmp$census_msoa))

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
png("/path/NDVI_locs_census_MSOA_20240904.png", width=1500, height=2000)

plot(msoa_merged, main="UK Biobank population with NDVI by MSOA", col=pal[msoa_merged$colours], lwd=0.5, cex.main=3)

legend(470000, 1000000, legend=c("0", "1 to 100", "101 to 200", "201 to 300", "301 to 400", "401 to 500", "501 to 600", "601 to 700", "700+"), col=pal, pch=15, cex=2.3, pt.cex=2.3, title="N partcipants per MSOA")

dev.off()

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

#R-squared for depression diagnosis
mod<-glm(phen$depression_diagnosis ~ phen$grs, family="binomial")
PseudoR2(mod, which="Nagelkerke")

######################################
#Phenotypic analyses
######################################

#Association between depression diagnosis and NDVI, adjusted for sex and age
mod<-lm(phen$NDVI_500m_mean ~ phen$depression_diagnosis + phen$sex + phen$age)
summary(mod)
nobs(mod)
confint.default(mod)

#Association between MHQ dep and NDVI, adjusted for sex and age
mod<-lm(phen$NDVI_500m_mean ~ phen$MHQ_dep + phen$sex + phen$age)
summary(mod)
nobs(mod)
confint.default(mod)

#Association between depression diagnosis and percentage greenspace, adjusted for sex and age
mod<-lm(phen$greenspace_percent_300m ~ phen$depression_diagnosis + phen$sex + phen$age)
summary(mod)
nobs(mod)
confint.default(mod)

#Association between MHQ dep and percentage greenspace, adjusted for sex and age
mod<-lm(phen$greenspace_percent_300m ~ phen$MHQ_dep + phen$sex + phen$age)
summary(mod)
nobs(mod)
confint.default(mod)

