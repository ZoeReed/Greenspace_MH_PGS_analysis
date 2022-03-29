######################################
##Script to conduct one sample MR for depression on NDVI
#The syntax was created by Zoe E Reed.
#The syntax was checked by Tim T. Morris
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################
library(AER)
library(stats)
library(simex)

######################################
#Read in data
######################################

#Phenotypic
load("/filepath/phenotypic_one_sample_MR.RData")

#PGS
dep_grs<-read.csv(file="/filepath/grs_dep_5e8.csv")

######################################
#Format data
######################################

#Merge datasets into a single data frame
all<-merge(phen, dep_grs, by.x="IID", by.y="id") 

all$sex<-as.factor(all$sex)
all$MHQ_dep<-as.factor(all$MHQ_dep)
all$depression_diagnosis<-as.factor(all$depression_diagnosis)
all$depression_non_cancer<-as.factor(all$depression_non_cancer)
all$depression_icd_primary<-as.factor(all$depression_icd_primary)
all$depression_icd_secondary<-as.factor(all$depression_icd_secondary)
all$Dep_MH_online<-as.factor(all$Dep_MH_online)

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

######################################
#Descriptives
######################################

str(all)

#NDVI
sum(!is.na(all$NDVI_500m_mean))
mean(all$NDVI_500m_mean, na.rm=T)
sd(all$NDVI_500m_mean, na.rm=T)

#Age
mean(all$age, na.rm=T)
sd(all$age, na.rm=T)
range(all$age, na.rm=T)

#Sex, 1 is male
table(all$sex)

#Depression diagnosis
table(all$depression_diagnosis)

#Depression MHQ
table(all$MHQ_dep)

#Depression MH online
table(all$Dep_MH_online)

#Depression ICD primary and secondary
table(all$depression_icd_primary)
table(all$depression_icd_secondary)
table(all$depression_icd_primary, all$depression_icd_secondary) 

#Depression non cancer codes
table(all$depression_non_cancer) #5.67%

######################################
#Standardise
######################################

#Standardise NDVI
all$NDVI_500m_mean<-scale(all$NDVI_500m_mean, center=T, scale=T)

#Standardise PGS
all$grs<-scale(all$grs, center=T, scale=T)

######################################
#Check associations
######################################

#NDVI and depression
summary(lm(all$NDVI_500m_mean ~ all$depression_diagnosis)) 
nobs(lm(all$NDVI_500m_mean ~ all$depression_diagnosis)) #238,306

summary(lm(all$NDVI_500m_mean ~ all$depression_diagnosis + all$age + all$sex))
nobs(lm(all$NDVI_500m_mean ~ all$depression_diagnosis + all$age + all$sex)) #238,306

#Sensitivity analysis for MHQ depression
summary(lm(all$NDVI_500m_mean ~ all$MHQ_dep)) 
nobs(lm(all$NDVI_500m_mean ~ all$MHQ_dep)) #64,740

summary(lm(all$NDVI_500m_mean ~ all$MHQ_dep + all$age + all$sex))
nobs(lm(all$NDVI_500m_mean ~ all$MHQ_dep + all$age + all$sex)) #64,740

#Depression and confounders (age, sex)
#Age (dep diagnois)
model1<-lm(all$age~all$depression_diagnosis)
summary(model1)

#Sex (dep diagnois)
chi_mod<-chisq.test(all$depression_diagnosis,all$sex)
chi_mod$observed

#Age (dep MHQ)
model1<-lm(all$age~all$MHQ_dep)
summary(model1)

#Sex (dep MHQ)
chi_mod<-chisq.test(all$MHQ_dep,all$sex)
chi_mod$observed

#PGS and confounders (age, sex)
#Age 
model2<-lm(all$grs~all$age)
summary(model1)

#Sex
t.test(all$grs~all$sex) #no difference

######################################
#Conduct 2SLS regression
######################################

tsls<-ivreg(NDVI_500m_mean ~ depression_diagnosis + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 + PC24 + PC25 | grs + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 + PC24 + PC25, data=all)

summary(tsls, diagnostics=T)
nobs(tsls) #238,306
confint.default(tsls)

#Sensitivity analysis with MHQ dep
tsls<-ivreg(NDVI_500m_mean ~ MHQ_dep + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 + PC24 + PC25 | grs + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 + PC24 + PC25, data=all)

summary(tsls, diagnostics=T)
nobs(tsls, diagnostics=T) #64,740
confint.default(tsls) 
