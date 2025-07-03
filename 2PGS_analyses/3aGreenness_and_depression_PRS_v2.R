######################################
##Script to test the association between the Depression PRS (at multiple thresholds) and the NDVI variable
#The syntax was created by Zoe E Reed.
#The syntax was checked by Gareth J. Griffith
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries
######################################

library(data.table)
library(plyr)
library(dplyr)
library(lme4)

######################################
#Load data and format
######################################

#Load data
load("/path/greenspace_biobank_final_census_20240904.RData")

#Load PCs data and match
pcs<-read.table("/path/data.pca1-40.plink.txt")
colnames(pcs)<-c("IID", "FID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", "PC33", "PC34", "PC35", "PC36", "PC37", "PC38", "PC39", "PC40")
vars<-c("IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25")
pcs<-pcs[,vars]

common.samples<-intersect(as.character(data$IID), as.character(pcs$IID))
data<-data[match(common.samples, data$IID),]
pcs<-pcs[match(common.samples, pcs$IID),]
data2<-left_join(data, pcs, by="IID")

#Load PRS data
polygenic_risk_score_5e8<-read.csv("/path/grs_dep_5e8.csv")
polygenic_risk_score_000001<-read.csv("/path/grs_dep_000001.csv")
polygenic_risk_score_00001<-read.csv("/path/grs_dep_00001.csv")
polygenic_risk_score_0001<-read.csv("/path/grs_dep_0001.csv")
polygenic_risk_score_001<-read.csv("/path/grs_dep_001.csv")
polygenic_risk_score_01<-read.csv("/path/grs_dep_01.csv")
polygenic_risk_score_05<-read.csv("/path/grs_dep_05.csv")
polygenic_risk_score_1<-read.csv("/path/grs_dep_1.csv")
polygenic_risk_score_2<-read.csv("/path/grs_dep_2.csv")
polygenic_risk_score_3<-read.csv("/path/grs_dep_3.csv")
polygenic_risk_score_4<-read.csv("/path/grs_dep_4.csv")
polygenic_risk_score_5<-read.csv("/path/grs_dep_5.csv")

#Standardise scores
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
rec_exc<-read.table("/path/data.combined_recommended.qctools.txt", header=F)
data2<-data2[!(data2$IID %in% rec_exc$V1), ]
dim(data2)

high_rel<-read.csv("/path/data.highly_relateds.qctools.txt", header=F)
data2<-data2[!(data2$IID %in% high_rel$V1), ]
dim(data2)

min_rel<-read.csv("/path/data.minimal_relateds.qctools.txt", header=F)
data2<-data2[!(data2$IID %in% min_rel$V1), ]
dim(data2)

non_white_brit<-read.csv("/path/data.non_white_british.qctools.txt", header=F)
data2<-data2[!(data2$IID %in% non_white_brit$V1), ]
dim(data2)

data2<-data2[!(is.na(data2$NDVI_500m_mean)), ]
dim(data2)

#Standardise variable
data2$NDVI_500m_mean<-scale(data2$NDVI_500m_mean, center=T, scale=T)

#Match to phenotypic dataset
common.samples<-intersect(as.character(data2$IID), as.character(polygenic_risk_score_5$id))
data2<-data2[match(common.samples, data2$IID),]

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
#Analyses
######################################

Greenness_Dep_PRS_lm_results<-list()
i<-0

for(risk_score in scores){
	i<-i+1
	tmp_data<-merge(data2, risk_score[,c("id", "z_risk_score")], by.x="IID", by.y="id")
	greenness_agg<-tmp_data %>% group_by(census_msoa) %>% 
	summarise_at(vars(sex, age, NDVI_500m_mean, z_risk_score, PC1, PC2, PC3, PC4, PC5,PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20, PC21, PC22, PC23, PC24, PC25), list(mean=mean))
	greenness_edited <- left_join(tmp_data, greenness_agg[c("census_msoa","z_risk_score_mean")])
	## List variables for model specs
	outcome <- "NDVI_500m_mean"
	exposure <- "z_risk_score"
	mundlak_exp <- "z_risk_score_mean"
	confs <- c("sex", "age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25")
	re <- "(1|census_msoa)"
	summ_lm <- glm(as.formula(paste(outcome,
									paste(c(exposure,
                                    		confs),
                                    	collapse = " + "),
                                	sep="~")),
               		data=greenness_edited)
    modelb<-lm(NDVI_500m_mean ~ z_risk_score, data=greenness_edited)
	rsq<-summary(modelb)$r.squared
	new<-as.data.frame(rsq)
	new$score<-names(scores[i])
	new$beta<-summary(summ_lm)$coefficients[2,1]
	new$p<-summary(summ_lm)$coefficients[2,4]
	new$Lower_CI<-confint.default(summ_lm)[2,1]
	new$Upper_CI<-confint.default(summ_lm)[2,2]
	new$N<-nobs(summ_lm)
	Greenness_Dep_PRS_lm_results<-rbind(Greenness_Dep_PRS_lm_results, new)
	}

write.csv(Greenness_Dep_PRS_lm_results, file="/path/Greenness_Dep_PRS_results_lm.csv", row.names=F, quote=F)

#MSOA with Mundlak formula
rm(tmp_data)
rm(greenness_edited)
rm(greenness_agg)
rm(outcome)
rm(exposure)
rm(mundlak_exp)
rm(confs)

Greenness_Dep_PRS_msoa_results_mundlak_risk_score<-list()
Greenness_Dep_PRS_msoa_results_mundlak_risk_score_mean<-list()
Residuals<-list()
i<-0

for(risk_score in scores){
	i<-i+1
	tmp_data<-merge(data2, risk_score[,c("id", "z_risk_score")], by.x="IID", by.y="id")
	greenness_agg<-tmp_data %>% group_by(census_msoa) %>% 
	summarise_at(vars(sex, age, NDVI_500m_mean, z_risk_score, PC1, PC2, PC3, PC4, PC5,PC6, PC7, PC8, PC9, PC10, PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20, PC21, PC22, PC23, PC24, PC25), list(mean=mean))
	greenness_edited <- left_join(tmp_data, greenness_agg[c("census_msoa","z_risk_score_mean")])
	#List variables for model specs
	outcome <- "NDVI_500m_mean"
	exposure <- "z_risk_score"
	mundlak_exp <- "z_risk_score_mean"
	confs <- c("sex", "age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25")
	re <- "(1|census_msoa)"
	summ_lm <- lmer(as.formula(paste(outcome,
									paste(c(exposure,
											mundlak_exp,
                                    		confs,
                                    		re),
                                    	collapse = " + "),
                                	sep="~")),
               		data=greenness_edited)
    modelb<-lm(NDVI_500m_mean ~ z_risk_score, data=greenness_edited)
    u0 <- ranef(summ_lm, condVar = TRUE) #extract level 2 residuals
	N<-nobs(summ_lm)
	new<-as.data.frame(N)
	new$score<-names(scores[i])
	new$beta<-summary(summ_lm)$coefficients[2,1]
	new$t<-summary(summ_lm)$coefficients[2,3]
	new$SE<-summary(summ_lm)$coefficients[2,2]
	new$Lower_CI<-new$beta-(1.96*new$SE)
	new$Upper_CI<-new$beta+(1.96*new$SE)
	coefs<-data.frame(coef(summary(summ_lm)))
	new$p<-2*(1-pnorm(abs(coefs$t.value[2])))
	new2<-as.data.frame(N)
	new2$score<-names(scores[i])
	new2$beta<-summary(summ_lm)$coefficients[3,1]
	new2$t<-summary(summ_lm)$coefficients[3,3]
	new2$SE<-summary(summ_lm)$coefficients[3,2]
	new2$Lower_CI<-new2$beta-(1.96*new2$SE)
	new2$Upper_CI<-new2$beta+(1.96*new2$SE)
	coefs2<-data.frame(coef(summary(summ_lm)))
	new2$p<-2*(1-pnorm(abs(coefs2$t.value[3])))
	Greenness_Dep_PRS_msoa_results_mundlak_risk_score<-rbind(Greenness_Dep_PRS_msoa_results_mundlak_risk_score, new)
	Greenness_Dep_PRS_msoa_results_mundlak_risk_score_mean<-rbind(Greenness_Dep_PRS_msoa_results_mundlak_risk_score_mean, new2)
	Residuals<-rbind(Residuals, u0)
}

write.csv(Greenness_Dep_PRS_msoa_results_mundlak_risk_score, file="/path/Greenness_Dep_PRS_MSOA_results_Mundlak_risk_score.csv", row.names=F, quote=F)

write.csv(Greenness_Dep_PRS_msoa_results_mundlak_risk_score_mean, file="/path/Greenness_Dep_PRS_MSOA_results_Mundlak_risk_score_mean.csv", row.names=F, quote=F)

