######################################
##Script to test association between schizophrenia PGS (at multiple thresholds) and percentage greenspace
#The syntax was created by Zoe E Reed.
#The syntax was checked by Gareth J. Griffith
#R Version 3.5.1
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################
library(data.table)
library(plyr)
library(dplyr)
library(readstata13)
library(lme4)

######################################
#Load data
######################################

#Percentage greenspace dataset
load("/filepath/greenspace_biobank_final_census_20200805.RData")

#Principal components
pcs<-read.table("/filepath/data.pca1-40.plink.txt")
colnames(pcs)<-c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", "PC33", "PC34", "PC35", "PC36", "PC37", "PC38", "PC39", "PC40")
vars<-c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25")
pcs<-pcs[,vars]

#Match datasets
common.samples<-intersect(as.character(green_data$IID), as.character(pcs$IID))
green_data<-green_data[match(common.samples, green_data$IID),]
pcs<-pcs[match(common.samples, pcs$IID),]
green_data2<-left_join(green_data, pcs, by="IID")

#Load PGS data
polygenic_risk_score_5e8<-read.csv("/filepath/grs_scz_5e8.csv")
polygenic_risk_score_000001<-read.csv("/filepath/grs_scz_000001.csv")
polygenic_risk_score_00001<-read.csv("/filepath/grs_scz_00001.csv")
polygenic_risk_score_0001<-read.csv("/filepath/grs_scz_0001.csv")
polygenic_risk_score_001<-read.csv("/filepath/grs_scz_001.csv")
polygenic_risk_score_01<-read.csv("/filepath/grs_scz_01.csv")
polygenic_risk_score_05<-read.csv("/filepath/grs_scz_05.csv")
polygenic_risk_score_1<-read.csv("/filepath/grs_scz_1.csv")
polygenic_risk_score_2<-read.csv("/filepath/grs_scz_2.csv")
polygenic_risk_score_3<-read.csv("/filepath/grs_scz_3.csv")
polygenic_risk_score_4<-read.csv("/filepath/grs_scz_4.csv")
polygenic_risk_score_5<-read.csv("/filepath/grs_scz_5.csv")

######################################
#Format data
######################################

#Standardise PGS
polygenic_risk_score_5e8$z_risk_score<-scale(polygenic_risk_score_5e8$grs, center=T, scale=T)
polygenic_risk_score_000001$z_risk_score<-scale(polygenic_risk_score_000001$grs, center=T, scale=T)
polygenic_risk_score_00001$z_risk_score<-scale(polygenic_risk_score_00001$grs, center=T, scale=T)
polygenic_risk_score_0001$z_risk_score<-scale(polygenic_risk_score_0001$grs, center=T, scale=T)
polygenic_risk_score_001$z_risk_score<-scale(polygenic_risk_score_001$grs, center=T, scale=T)
polygenic_risk_score_01$z_risk_score<-scale(polygenic_risk_score_01$grs, center=T, scale=T)
polygenic_risk_score_05$z_risk_score<-scale(polygenic_risk_score_05$grs, center=T, scale=T)
polygenic_risk_score_1$z_risk_score<-scale(polygenic_risk_score_1$grs, center=T, scale=T)
polygenic_risk_score_2$z_risk_score<-scale(polygenic_risk_score_2$grs, center=T, scale=T)
polygenic_risk_score_3$z_risk_score<-scale(polygenic_risk_score_3$grs, center=T, scale=T)
polygenic_risk_score_4$z_risk_score<-scale(polygenic_risk_score_4$grs, center=T, scale=T)
polygenic_risk_score_5$z_risk_score<-scale(polygenic_risk_score_5$grs, center=T, scale=T)

#Remove exclusions for genetic data
rec_exc<-read.table("/filepath/data.combined_recommended.qctools.txt", header=F)
green_data2<-green_data2[!(green_data2$IID %in% rec_exc$V1), ]
dim(green_data2)

high_rel<-read.csv("/filepath/data.highly_relateds.qctools.txt", header=F)
green_data2<-green_data2[!(green_data2$IID %in% high_rel$V1), ]
dim(green_data2)

min_rel<-read.csv("/filepath/data.minimal_relateds.qctools.txt", header=F)
green_data2<-green_data2[!(green_data2$IID %in% min_rel$V1), ]
dim(green_data2)

non_white_brit<-read.csv("/filepath/data.non_white_british.qctools.txt", header=F)
green_data2<-green_data2[!(green_data2$IID %in% non_white_brit$V1), ]
dim(green_data2)

green_data2<-green_data2[!(is.na(green_data2$greenspace_percent_300m)), ]
dim(green_data2) #N=293,922 participants

#Standardise variable
green_data2$greenspace_percent_300m<-scale(green_data2$greenspace_percent_300m, center=T, scale=T)

#Match PGS to phenotypic dataset
common.samples<-intersect(as.character(green_data2$IID), as.character(polygenic_risk_score_5$id))
green_data2<-green_data2[match(common.samples, green_data2$IID),]

polygenic_risk_score_5e8<-polygenic_risk_score_5e8[match(common.samples, polygenic_risk_score_5e8$id),]
polygenic_risk_score_000001<-polygenic_risk_score_000001[match(common.samples, polygenic_risk_score_000001$id),]
polygenic_risk_score_00001<-polygenic_risk_score_00001[match(common.samples, polygenic_risk_score_00001$id),]
polygenic_risk_score_0001<-polygenic_risk_score_0001[match(common.samples, polygenic_risk_score_0001$id),]
polygenic_risk_score_001<-polygenic_risk_score_001[match(common.samples, polygenic_risk_score_001$id),]
polygenic_risk_score_01<-polygenic_risk_score_01[match(common.samples, polygenic_risk_score_01$id),]
polygenic_risk_score_05<-polygenic_risk_score_05[match(common.samples, polygenic_risk_score_05$id),]
polygenic_risk_score_1<-polygenic_risk_score_1[match(common.samples, polygenic_risk_score_1$id),]
polygenic_risk_score_2<-polygenic_risk_score_2[match(common.samples, polygenic_risk_score_2$id),]
polygenic_risk_score_3<-polygenic_risk_score_3[match(common.samples, polygenic_risk_score_3$id),]
polygenic_risk_score_4<-polygenic_risk_score_4[match(common.samples, polygenic_risk_score_4$id),]
polygenic_risk_score_5<-polygenic_risk_score_5[match(common.samples, polygenic_risk_score_5$id),]

scores<-list(polygenic_risk_score_5e8, polygenic_risk_score_000001, polygenic_risk_score_00001, polygenic_risk_score_0001, polygenic_risk_score_001, polygenic_risk_score_01, polygenic_risk_score_05, polygenic_risk_score_1, polygenic_risk_score_2, polygenic_risk_score_3, polygenic_risk_score_4, polygenic_risk_score_5)

names(scores)<-c("polygenic_score_5e8", "polygenic_score_000001", "polygenic_score_00001", "polygenic_score_0001", "polygenic_score_001", "polygenic_score_01", "polygenic_score_05", "polygenic_score_1", "polygenic_score_2", "polygenic_score_3", "polygenic_score_4", "polygenic_score_5")

######################################
#Analysis
######################################

##Greenspace~Scz PRS
Greenspace_Scz_PRS_results<-list()
i<-0

for(risk_score in scores){
	i<-i+1
	model<-lm(green_data2$greenspace_percent_300m ~ risk_score$z_risk_score + green_data2$sex + green_data2$age + green_data2$PC1 + green_data2$PC2 + green_data2$PC3 + green_data2$PC4 + green_data2$PC5 + green_data2$PC6 + green_data2$PC7 + green_data2$PC8 + green_data2$PC9 + green_data2$PC10 + green_data2$PC11 + green_data2$PC12 + green_data2$PC13 + green_data2$PC14 + green_data2$PC15 + green_data2$PC16 + green_data2$PC17 + green_data2$PC18 + green_data2$PC19 + green_data2$PC20 + green_data2$PC21 + green_data2$PC22 + green_data2$PC23 + green_data2$PC24 + green_data2$PC25)
	modelb<-lm(green_data2$greenspace_percent_300m ~ risk_score$z_risk_score)
	rsq<-summary(modelb)$r.squared
	new<-as.data.frame(rsq)
	new$score<-names(scores[i])
	new$beta<-summary(model)$coefficients[2,1]
	new$p<-summary(model)$coefficients[2,4]
	new$Lower_CI<-confint.default(model)[2,1]
	new$Upper_CI<-confint.default(model)[2,2]
	new$N<-nobs(model)
	Greenspace_Scz_PRS_results<-rbind(Greenspace_Scz_PRS_results, new)
}

write.table(Greenspace_Scz_PRS_results, file="/filepath/Greenspace_Scz_PRS_results.csv", row.names=F, append=F, col.names=T, sep="\t", quote=F)

#MSOA with random intercept
Greenspace_Scz_PRS_msoa_results<-list()
i<-0

for(risk_score in scores){
	i<-i+1
	model<-lmer(green_data2$greenspace_percent_300m ~ risk_score$z_risk_score + green_data2$sex + green_data2$age + green_data2$PC1 + green_data2$PC2 + green_data2$PC3 + green_data2$PC4 + green_data2$PC5 + green_data2$PC6 + green_data2$PC7 + green_data2$PC8 + green_data2$PC9 + green_data2$PC10 + green_data2$PC11 + green_data2$PC12 + green_data2$PC13 + green_data2$PC14 + green_data2$PC15 + green_data2$PC16 + green_data2$PC17 + green_data2$PC18 + green_data2$PC19 + green_data2$PC20 + green_data2$PC21 + green_data2$PC22 + green_data2$PC23 + green_data2$PC24 + green_data2$PC25 + (1|green_data2$census_msoa))
	N<-nobs(model)
	new<-as.data.frame(N)
	new$score<-names(scores[i])
	new$beta<-summary(model)$coefficients[2,1]
	new$t<-summary(model)$coefficients[2,3]
	new$SE<-summary(model)$coefficients[2,2]
	new$Lower_CI<-new$beta-(1.96*new$SE)
	new$Upper_CI<-new$beta+(1.96*new$SE)
	coefs<-data.frame(coef(summary(model)))
	new$p<-2*(1-pnorm(abs(coefs$t.value[2])))
	Greenspace_Scz_PRS_msoa_results<-rbind(Greenspace_Scz_PRS_msoa_results, new)
}

write.table(Greenspace_Scz_PRS_msoa_results, file="/filepath/Greenspace_Scz_PRS_MSOA_results.csv", row.names=F, append=F, col.names=T, sep="\t", quote=F)

#VCP for MSOA
Greenspace_Scz_PRS_msoa_vpc<-list()
i<-0

for(risk_score in scores){
	i<-i+1
	model<-lmer(green_data2$greenspace_percent_300m ~ risk_score$z_risk_score + green_data2$sex + green_data2$age + green_data2$PC1 + green_data2$PC2 + green_data2$PC3 + green_data2$PC4 + green_data2$PC5 + green_data2$PC6 + green_data2$PC7 + green_data2$PC8 + green_data2$PC9 + green_data2$PC10 + green_data2$PC11 + green_data2$PC12 + green_data2$PC13 + green_data2$PC14 + green_data2$PC15 + green_data2$PC16 + green_data2$PC17 + green_data2$PC18 + green_data2$PC19 + green_data2$PC20 + green_data2$PC21 + green_data2$PC22 + green_data2$PC23 + green_data2$PC24 + green_data2$PC25 + (1|green_data2$census_msoa))
	vpc<-VarCorr(model) %>%
			as.data.frame() %>%
			mutate(icc=vcov/sum(vcov)) %>%
			select(grp, icc)
	new<-as.data.frame(vpc$icc[1])
	new$score<-names(scores[i])
	Greenspace_Scz_PRS_msoa_vpc<-rbind(Greenspace_Scz_PRS_msoa_vpc, new)
}

write.table(Greenspace_Scz_PRS_msoa_vpc, file="/filepath/Greenspace_Scz_PRS_MSOA_VPC.csv", row.names=F, append=F, col.names=T, sep="\t", quote=F)

#LAD with random intercept
Greenspace_Scz_PRS_lad_results<-list()
i<-0

for(risk_score in scores){
	i<-i+1
	model<-lmer(green_data2$greenspace_percent_300m ~ risk_score$z_risk_score + green_data2$sex + green_data2$age + green_data2$PC1 + green_data2$PC2 + green_data2$PC3 + green_data2$PC4 + green_data2$PC5 + green_data2$PC6 + green_data2$PC7 + green_data2$PC8 + green_data2$PC9 + green_data2$PC10 + green_data2$PC11 + green_data2$PC12 + green_data2$PC13 + green_data2$PC14 + green_data2$PC15 + green_data2$PC16 + green_data2$PC17 + green_data2$PC18 + green_data2$PC19 + green_data2$PC20 + green_data2$PC21 + green_data2$PC22 + green_data2$PC23 + green_data2$PC24 + green_data2$PC25 + (1|green_data2$census_lad))
	N<-nobs(model)
	new<-as.data.frame(N)
	new$score<-names(scores[i])
	new$beta<-summary(model)$coefficients[2,1]
	new$t<-summary(model)$coefficients[2,3]
	new$SE<-summary(model)$coefficients[2,2]
	new$Lower_CI<-new$beta-(1.96*new$SE)
	new$Upper_CI<-new$beta+(1.96*new$SE)
	coefs<-data.frame(coef(summary(model)))
	new$p<-2*(1-pnorm(abs(coefs$t.value[2])))
	Greenspace_Scz_PRS_lad_results<-rbind(Greenspace_Scz_PRS_lad_results, new)
}

write.table(Greenspace_Scz_PRS_lad_results, file="/filepath/Greenspace_Scz_PRS_LAD_results.csv", row.names=F, append=F, col.names=T, sep="\t", quote=F)

#VCP for LAD
Greenspace_Scz_PRS_lad_vpc<-list()
i<-0

for(risk_score in scores){
	i<-i+1
	model<-lmer(green_data2$greenspace_percent_300m ~ risk_score$z_risk_score + green_data2$sex + green_data2$age + green_data2$PC1 + green_data2$PC2 + green_data2$PC3 + green_data2$PC4 + green_data2$PC5 + green_data2$PC6 + green_data2$PC7 + green_data2$PC8 + green_data2$PC9 + green_data2$PC10 + green_data2$PC11 + green_data2$PC12 + green_data2$PC13 + green_data2$PC14 + green_data2$PC15 + green_data2$PC16 + green_data2$PC17 + green_data2$PC18 + green_data2$PC19 + green_data2$PC20 + green_data2$PC21 + green_data2$PC22 + green_data2$PC23 + green_data2$PC24 + green_data2$PC25 + (1|green_data2$census_lad))
	vpc<-VarCorr(model) %>%
			as.data.frame() %>%
			mutate(icc=vcov/sum(vcov)) %>%
			select(grp, icc)
	new<-as.data.frame(vpc$icc[1])
	new$score<-names(scores[i])
	Greenspace_Scz_PRS_lad_vpc<-rbind(Greenspace_Scz_PRS_lad_vpc, new)
}

write.table(Greenspace_Scz_PRS_lad_vpc, file="/filepath/Greenspace_Scz_PRS_LAD_VPC.csv", row.names=F, append=F, col.names=T, sep="\t", quote=F)
