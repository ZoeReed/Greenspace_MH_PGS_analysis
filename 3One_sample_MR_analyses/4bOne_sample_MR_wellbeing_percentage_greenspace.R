######################################
##Script to conduct one sample MR for wellbeing on percentage greenspace
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
well_grs<-read.csv(file="/filepath/grs_well_5e8.csv")

well_grs_1x10min05<-read.csv(file="/filepath/grs_well_00001.csv")

colnames(well_grs_1x10min05)<-c("id", "grs_1x10min05")

######################################
#Format data
######################################

#Merge datasets into a single data frame
tmp<-merge(phen, well_grs, by.x="IID", by.y="id") 

all<-merge(tmp, well_grs_1x10min05, by.x="IID", by.y="id") 


all$sex<-as.factor(all$sex)

#Treat wellbeing as continuous phenotype
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

######################################
#Descriptives
######################################

str(all)

#Percentage greenspace
sum(!is.na(all$greenspace_percent_300m))
mean(all$greenspace_percent_300m, na.rm=T)
sd(all$greenspace_percent_300m, na.rm=T)

#Wellbeing
sum(!is.na(all$wellbeing))
mean(all$wellbeing, na.rm=T)
sd(all$wellbeing, na.rm=T)

#Age
mean(all$age, na.rm=T)
sd(all$age, na.rm=T)
range(all$age, na.rm=T)

#Sex, 1 is male
table(all$sex)

######################################
#Standardise
######################################

#Standardise percentage greenspace
all$greenspace_percent_300m<-scale(all$greenspace_percent_300m, center=T, scale=T)

#Standardise PGS
all$grs<-scale(all$grs, center=T, scale=T)

######################################
#Check associations
######################################

#Percentage greenspace and depression
summary(lm(all$greenspace_percent_300m ~ all$wellbeing)) 
nobs(lm(all$greenspace_percent_300m ~ all$wellbeing)) #97,099

summary(lm(all$greenspace_percent_300m ~ all$wellbeing + all$age + all$sex))
nobs(lm(all$greenspace_percent_300m ~ all$wellbeing + all$age + all$sex)) #97,099

#Wellbeing and confounders (age, sex)
#Age
model1<-lm(all$wellbeing~all$age)
summary(model1)

#Sex
t.test(all$wellbeing~all$sex)

#PGS and confounders (age, sex)
#Age 
model2<-lm(all$grs~all$age)
summary(model1)

#Sex
t.test(all$grs~all$sex)

#Age (other PGS)
model2<-lm(all$grs_1x10min05~all$age)
summary(model1)

#Sex (other PGS)
t.test(all$grs_1x10min05~all$sex)

######################################
#Conduct 2SLS regression
######################################

tsls<-ivreg(greenspace_percent_300m ~ wellbeing + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 + PC24 + PC25 | grs + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 + PC24 + PC25, data=all)

summary(tsls, diagnostics=T)
nobs(tsls) #97,099
confint.default(tsls)

#Sensitivity analysis with grs_1x10min05
tsls<-ivreg(greenspace_percent_300m ~ wellbeing + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 + PC24 + PC25 | grs_1x10min05 + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21 + PC22 + PC23 + PC24 + PC25, data=all)

summary(tsls, diagnostics=T)
nobs(tsls, diagnostics=T) #97,099
confint.default(tsls)
