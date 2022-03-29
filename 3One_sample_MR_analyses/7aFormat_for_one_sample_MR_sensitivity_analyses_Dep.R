######################################
##Script to format data for one sample MR sensitivity analyses for depression
#The syntax was created by Zoe E Reed.
#The syntax was checked by ....
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################
library(dplyr)
library(data.table)
library(stringr)

######################################
#Read in data
######################################

#Genotype data
Dep_gen<-read.table("/filepath/Dep_snps.raw", header=T)

#Phenotypic
load("/filepath/phenotypic_one_sample_MR.RData")

######################################
#Format data
######################################

phen$sex<-as.factor(phen$sex)

######################################
#Load PCs
######################################
pcs_all<-read.table("/filepath/data.pca1-40.plink.txt")
colnames(pcs_all)<-c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", "PC33", "PC34", "PC35", "PC36", "PC37", "PC38", "PC39", "PC40")
vars<-c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25")
pcs<-pcs_all[,vars]

#Merge data
tmp<-phen
phen<-merge(tmp, pcs, by="IID", all.x=T)
all_data<-merge(phen, Dep_gen, by="IID", all.x=T)

######################################
#Exclusions
######################################

#Make sure all withdrawn are removed
withdrawn<-read.csv("/filepath/w21829_20200204.csv", header=F)
all_data<-all_data[!(all_data$eid %in% withdrawn$V1), ]
dim(all_data)

rec_exc<-read.table("/filepath/data.combined_recommended.qctools.txt", header=F)
all_data<-all_data[!(all_data$IID %in% rec_exc$V1), ]
dim(all_data)

high_rel<-read.csv("/filepath/data.highly_relateds.qctools.txt", header=F)
all_data<-all_data[!(all_data$IID %in% high_rel$V1), ]
dim(all_data)

min_rel<-read.csv("/filepath/data.minimal_relateds.qctools.txt", header=F)
all_data<-all_data[!(all_data$IID %in% min_rel$V1), ]
dim(all_data)

non_white_brit<-read.csv("/filepath/data.non_white_british.qctools.txt", header=F)
all_data<-all_data[!(all_data$IID %in% non_white_brit$V1), ]
dim(all_data)

#Restrict to only those with genetic data (use one SNP as example)
all_data<-all_data[!is.na(all_data$rs5758265_A),]
dim(all_data) #336,997

######################################
#Standardise
######################################

#Standardise NDVI
all_data$NDVI_500m_mean<-scale(all_data$NDVI_500m_mean, center=T, scale=T)

#Standardise percentage greenspace
all_data$greenspace_percent_300m<-scale(all_data$greenspace_percent_300m, center=T, scale=T)

########################
#Format data for one sample MR sensitivity analyses
########################

#NDVI dataset (Dep diag and MHQ as exposures)
#Summary statistics of genetic associations of Depression and outcomes
#Will adjust for sex, age, 25 genetic principal components

n <- nrow(all_data)
size <- ncol(all_data[,-(1:50)]) #2nd number is number of columns up to SNP data
G.1 <- as.matrix(all_data[,-(1:50)]) #2nd number is number of columns up to SNP data
X.1 <- as.matrix(as.numeric(all_data[,8])) #Exposure X.1 is dep diag 
X.2 <- as.matrix(as.numeric(all_data[,7])) #Exposure X.2 is MHQ dep
Y <- as.matrix(all_data[,19]) #Outcome 

#Covariates: sex, age, 25 genetic principal components
sex<-as.matrix(all_data[,5])
age<-as.matrix(all_data[,4])

pc1 <- as.matrix(all_data[,21])     
pc2 <- as.matrix(all_data[,22]) 
pc3 <- as.matrix(all_data[,23]) 
pc4 <- as.matrix(all_data[,24]) 
pc5 <- as.matrix(all_data[,25]) 
pc6 <- as.matrix(all_data[,26]) 
pc7 <- as.matrix(all_data[,27]) 
pc8 <- as.matrix(all_data[,28]) 
pc9 <- as.matrix(all_data[,29]) 
pc10 <- as.matrix(all_data[,30]) 
pc11 <- as.matrix(all_data[,31])     
pc12 <- as.matrix(all_data[,32]) 
pc13 <- as.matrix(all_data[,33]) 
pc14 <- as.matrix(all_data[,34]) 
pc15 <- as.matrix(all_data[,35]) 
pc16 <- as.matrix(all_data[,36]) 
pc17 <- as.matrix(all_data[,37]) 
pc18 <- as.matrix(all_data[,38]) 
pc19 <- as.matrix(all_data[,39]) 
pc20 <- as.matrix(all_data[,40]) 
pc21 <- as.matrix(all_data[,41])     
pc22 <- as.matrix(all_data[,42]) 
pc23 <- as.matrix(all_data[,43]) 
pc24 <- as.matrix(all_data[,44]) 
pc25 <- as.matrix(all_data[,45]) 

#Depression Diagnosis
#X~G
XGdata1        = data.frame(X.1, G.1,
                            sex, age, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25)   

FIT2           = summary(lm(XGdata1))                      
BetaXG         = FIT2$coef[-1,1]                              
BetaXG         = head(BetaXG, n= -27)                           
seBetaXG       = FIT2$coef[-1,2]                                
seBetaXG       = head(seBetaXG, n= -27)                           
Fbar           = mean((BetaXG^2)/(seBetaXG^2))                 

#Y~G  
YGdata        = data.frame(Y, G.1,
                           sex, age, 
                           pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                           pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                           pc21, pc22, pc23, pc24, pc25)   

FIT3           = summary(lm(YGdata))                            
BetaYG         = FIT3$coef[-1,1]                              
BetaYG         = head(BetaYG, n= -27)                           
seBetaYG       = FIT3$coef[-1,2]                              
seBetaYG       = head(seBetaYG, n= -27)                           

#Y~G+X 
YXGdata        = data.frame(Y, X.1, G.1,
                            sex, age,
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25) 

FIT4             = summary(lm(YXGdata))                       
alphahatstar     = FIT4$coef[-c(1,2),1]                       
alphahatstar     = head(alphahatstar, n= -27)                          
se.alphahatstar  = FIT4$coef[-c(1,2),2]
se.alphahatstar  = head(se.alphahatstar, n= -27)                           
betastar         = FIT4$coef[2,1] 
se.betastar  	 = FIT4$coef[2,2] 

#Summary statistics
# X - G
BetaXG_data<-as.matrix(BetaXG)
colnames(BetaXG_data) <- c("BetaXG")
seBetaXG_data<-as.matrix(seBetaXG)
colnames(seBetaXG_data) <- c("seBetaXG")

#Y - G
BetaYG_data<-as.matrix(BetaYG)
colnames(BetaYG_data) <- c("BetaYG")
seBetaYG_data<-as.matrix(seBetaYG)
colnames(seBetaYG_data) <- c("seBetaYG")

#Y - X - G
alphahatstar_data<-as.matrix(alphahatstar)
colnames(alphahatstar_data) <- c("alphahatstar")
se.alphahatstar_data<-as.matrix(se.alphahatstar)
colnames(se.alphahatstar_data) <- c("se.alphahatstar")

betastar_data<-as.matrix(betastar)
colnames(betastar_data) <- c("betastar")
se.betastar_data<-as.matrix(se.betastar)
colnames(se.betastar_data) <- c("se.betastar")

#Combine data
data<-cbind(BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)

SNP<-as.matrix(row.names(data))
colnames(SNP) <- c("SNP_ukb")

data<-cbind(SNP, BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)
rownames(data) <- c()

########################
#Save data and remove for next part
########################

write.csv(data, "/filepath/all_data_dep_diag_adjusted_summary_data_ukb.csv", row.names = F)

rm(XGdata1)
rm(FIT2)
rm(BetaXG)
rm(seBetaXG)
rm(Fbar)
rm(YGdata)
rm(FIT3)
rm(BetaYG)
rm(seBetaYG)
rm(YXGdata)
rm(FIT4)
rm(alphahatstar)
rm(se.alphahatstar)
rm(betastar)
rm(se.betastar)
rm(BetaXG_data)
rm(seBetaXG_data)
rm(BetaYG_data)
rm(seBetaYG_data)
rm(alphahatstar_data)
rm(se.alphahatstar_data)
rm(betastar_data)
rm(se.betastar_data)
rm(data)
rm(SNP)

#Depression MHQ
# X~G
XGdata1        = data.frame(X.2, G.1,
                            sex, age, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25)   

FIT2           = summary(lm(XGdata1))                      
BetaXG         = FIT2$coef[-1,1]                              
BetaXG         = head(BetaXG, n= -27)                           
seBetaXG       = FIT2$coef[-1,2]                                
seBetaXG       = head(seBetaXG, n= -27)                           
Fbar           = mean((BetaXG^2)/(seBetaXG^2))                 

# Y~G  
YGdata        = data.frame(Y, G.1,
                           sex, age, 
                           pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                           pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                           pc21, pc22, pc23, pc24, pc25)   

FIT3           = summary(lm(YGdata))                            
BetaYG         = FIT3$coef[-1,1]                              
BetaYG         = head(BetaYG, n= -27)                           
seBetaYG       = FIT3$coef[-1,2]                              
seBetaYG       = head(seBetaYG, n= -27)                           

# Y~G+X 
YXGdata        = data.frame(Y, X.2, G.1,
                            sex, age,
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25) 


FIT4             = summary(lm(YXGdata))                       
alphahatstar     = FIT4$coef[-c(1,2),1]                       
alphahatstar     = head(alphahatstar, n= -27)                          
se.alphahatstar  = FIT4$coef[-c(1,2),2]
se.alphahatstar  = head(se.alphahatstar, n= -27)                           
betastar         = FIT4$coef[2,1] 
se.betastar  	 = FIT4$coef[2,2] 

#Summary statistics
#X - G
BetaXG_data<-as.matrix(BetaXG)
colnames(BetaXG_data) <- c("BetaXG")
seBetaXG_data<-as.matrix(seBetaXG)
colnames(seBetaXG_data) <- c("seBetaXG")

#Y - G
BetaYG_data<-as.matrix(BetaYG)
colnames(BetaYG_data) <- c("BetaYG")
seBetaYG_data<-as.matrix(seBetaYG)
colnames(seBetaYG_data) <- c("seBetaYG")

#Y - X - G
alphahatstar_data<-as.matrix(alphahatstar)
colnames(alphahatstar_data) <- c("alphahatstar")
se.alphahatstar_data<-as.matrix(se.alphahatstar)
colnames(se.alphahatstar_data) <- c("se.alphahatstar")

betastar_data<-as.matrix(betastar)
colnames(betastar_data) <- c("betastar")
se.betastar_data<-as.matrix(se.betastar)
colnames(se.betastar_data) <- c("se.betastar")

#Combine data
data<-cbind(BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)

SNP<-as.matrix(row.names(data))
colnames(SNP) <- c("SNP_ukb")

data<-cbind(SNP, BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)
rownames(data) <- c()

########################
#Save data and remove for next part
########################
write.csv(data, "/filepath/all_data_MHQ_dep_adjusted_summary_data_ukb.csv", row.names = F)

rm(XGdata1)
rm(FIT2)
rm(BetaXG)
rm(seBetaXG)
rm(Fbar)
rm(YGdata)
rm(FIT3)
rm(BetaYG)
rm(seBetaYG)
rm(YXGdata)
rm(FIT4)
rm(alphahatstar)
rm(se.alphahatstar)
rm(betastar)
rm(se.betastar)
rm(BetaXG_data)
rm(seBetaXG_data)
rm(BetaYG_data)
rm(seBetaYG_data)
rm(alphahatstar_data)
rm(se.alphahatstar_data)
rm(betastar_data)
rm(se.betastar_data)
rm(data)
rm(SNP)

#Percentage greenspace (Dep diag and MHQ as exposures)
#NDVI dataset (Dep diag and MHQ as exposures)
#Summary statistics of genetic associations of Depression and outcomes
#Will adjust for sex, age, 25 genetic principal components

n <- nrow(all_data)
size <- ncol(all_data[,-(1:50)]) #2nd number is number of columns up to SNP data
G.1 <- as.matrix(all_data[,-(1:50)]) #2nd number is number of columns up to SNP data
X.1 <- as.matrix(as.numeric(all_data[,8])) #Exposure X.1 is dep diag 
X.2 <- as.matrix(as.numeric(all_data[,7])) #Exposure X.2 is MHQ dep
Y <- as.matrix(all_data[,18]) #Outcome 

#Covariates: sex, age, 25 genetic principal components
sex<-as.matrix(all_data[,5])
age<-as.matrix(all_data[,4])

pc1 <- as.matrix(all_data[,21])     
pc2 <- as.matrix(all_data[,22]) 
pc3 <- as.matrix(all_data[,23]) 
pc4 <- as.matrix(all_data[,24]) 
pc5 <- as.matrix(all_data[,25]) 
pc6 <- as.matrix(all_data[,26]) 
pc7 <- as.matrix(all_data[,27]) 
pc8 <- as.matrix(all_data[,28]) 
pc9 <- as.matrix(all_data[,29]) 
pc10 <- as.matrix(all_data[,30]) 
pc11 <- as.matrix(all_data[,31])     
pc12 <- as.matrix(all_data[,32]) 
pc13 <- as.matrix(all_data[,33]) 
pc14 <- as.matrix(all_data[,34]) 
pc15 <- as.matrix(all_data[,35]) 
pc16 <- as.matrix(all_data[,36]) 
pc17 <- as.matrix(all_data[,37]) 
pc18 <- as.matrix(all_data[,38]) 
pc19 <- as.matrix(all_data[,39]) 
pc20 <- as.matrix(all_data[,40]) 
pc21 <- as.matrix(all_data[,41])     
pc22 <- as.matrix(all_data[,42]) 
pc23 <- as.matrix(all_data[,43]) 
pc24 <- as.matrix(all_data[,44]) 
pc25 <- as.matrix(all_data[,45]) 

#Depression Diagnosis
#X~G
XGdata1        = data.frame(X.1, G.1,
                            sex, age, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25)   

FIT2           = summary(lm(XGdata1))                      
BetaXG         = FIT2$coef[-1,1]                              
BetaXG         = head(BetaXG, n= -27)                           
seBetaXG       = FIT2$coef[-1,2]                                
seBetaXG       = head(seBetaXG, n= -27)                           
Fbar           = mean((BetaXG^2)/(seBetaXG^2))                 

#Y~G  
YGdata        = data.frame(Y, G.1,
                           sex, age, 
                           pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                           pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                           pc21, pc22, pc23, pc24, pc25)   

FIT3           = summary(lm(YGdata))                            
BetaYG         = FIT3$coef[-1,1]                              
BetaYG         = head(BetaYG, n= -27)                           
seBetaYG       = FIT3$coef[-1,2]                              
seBetaYG       = head(seBetaYG, n= -27)                           

#Y~G+X 
YXGdata        = data.frame(Y, X.1, G.1,
                            sex, age,
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25) 

FIT4             = summary(lm(YXGdata))                       
alphahatstar     = FIT4$coef[-c(1,2),1]                       
alphahatstar     = head(alphahatstar, n= -27)                          
se.alphahatstar  = FIT4$coef[-c(1,2),2]
se.alphahatstar  = head(se.alphahatstar, n= -27)                           
betastar         = FIT4$coef[2,1] 
se.betastar  	 = FIT4$coef[2,2] 

#Summary statistics
#X - G
BetaXG_data<-as.matrix(BetaXG)
colnames(BetaXG_data) <- c("BetaXG")
seBetaXG_data<-as.matrix(seBetaXG)
colnames(seBetaXG_data) <- c("seBetaXG")

#Y - G
BetaYG_data<-as.matrix(BetaYG)
colnames(BetaYG_data) <- c("BetaYG")
seBetaYG_data<-as.matrix(seBetaYG)
colnames(seBetaYG_data) <- c("seBetaYG")

#Y - X - G
alphahatstar_data<-as.matrix(alphahatstar)
colnames(alphahatstar_data) <- c("alphahatstar")
se.alphahatstar_data<-as.matrix(se.alphahatstar)
colnames(se.alphahatstar_data) <- c("se.alphahatstar")

betastar_data<-as.matrix(betastar)
colnames(betastar_data) <- c("betastar")
se.betastar_data<-as.matrix(se.betastar)
colnames(se.betastar_data) <- c("se.betastar")

#Combine data
data<-cbind(BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)

SNP<-as.matrix(row.names(data))
colnames(SNP) <- c("SNP_ukb")

data<-cbind(SNP, BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)
rownames(data) <- c()

########################
#Save data and remove for next part
########################

write.csv(data, "/filepath/Percentage_GS_data_dep_diag_adjusted_summary_data_ukb.csv", row.names = F)

rm(XGdata1)
rm(FIT2)
rm(BetaXG)
rm(seBetaXG)
rm(Fbar)
rm(YGdata)
rm(FIT3)
rm(BetaYG)
rm(seBetaYG)
rm(YXGdata)
rm(FIT4)
rm(alphahatstar)
rm(se.alphahatstar)
rm(betastar)
rm(se.betastar)
rm(BetaXG_data)
rm(seBetaXG_data)
rm(BetaYG_data)
rm(seBetaYG_data)
rm(alphahatstar_data)
rm(se.alphahatstar_data)
rm(betastar_data)
rm(se.betastar_data)
rm(data)
rm(SNP)

#Depression MHQ
#X~G
XGdata1        = data.frame(X.2, G.1,
                            sex, age, 
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25)   

FIT2           = summary(lm(XGdata1))                      
BetaXG         = FIT2$coef[-1,1]                              
BetaXG         = head(BetaXG, n= -27)                           
seBetaXG       = FIT2$coef[-1,2]                                
seBetaXG       = head(seBetaXG, n= -27)                           
Fbar           = mean((BetaXG^2)/(seBetaXG^2))                 

#Y~G  
YGdata        = data.frame(Y, G.1,
                           sex, age, 
                           pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                           pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                           pc21, pc22, pc23, pc24, pc25)   

FIT3           = summary(lm(YGdata))                            
BetaYG         = FIT3$coef[-1,1]                              
BetaYG         = head(BetaYG, n= -27)                           
seBetaYG       = FIT3$coef[-1,2]                              
seBetaYG       = head(seBetaYG, n= -27)                           

#Y~G+X 
YXGdata        = data.frame(Y, X.2, G.1,
                            sex, age,
                            pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10,
                            pc11, pc12, pc13, pc14, pc15, pc16, pc17, pc18, pc19, pc20,
                            pc21, pc22, pc23, pc24, pc25) 

FIT4             = summary(lm(YXGdata))                       
alphahatstar     = FIT4$coef[-c(1,2),1]                       
alphahatstar     = head(alphahatstar, n= -27)                          
se.alphahatstar  = FIT4$coef[-c(1,2),2]
se.alphahatstar  = head(se.alphahatstar, n= -27)                           
betastar         = FIT4$coef[2,1] 
se.betastar  	 = FIT4$coef[2,2] 

#Summary statistics
#X- G
BetaXG_data<-as.matrix(BetaXG)
colnames(BetaXG_data) <- c("BetaXG")
seBetaXG_data<-as.matrix(seBetaXG)
colnames(seBetaXG_data) <- c("seBetaXG")

#Y - G
BetaYG_data<-as.matrix(BetaYG)
colnames(BetaYG_data) <- c("BetaYG")
seBetaYG_data<-as.matrix(seBetaYG)
colnames(seBetaYG_data) <- c("seBetaYG")

#Y - X - G
alphahatstar_data<-as.matrix(alphahatstar)
colnames(alphahatstar_data) <- c("alphahatstar")
se.alphahatstar_data<-as.matrix(se.alphahatstar)
colnames(se.alphahatstar_data) <- c("se.alphahatstar")

betastar_data<-as.matrix(betastar)
colnames(betastar_data) <- c("betastar")
se.betastar_data<-as.matrix(se.betastar)
colnames(se.betastar_data) <- c("se.betastar")

#Combine data
data<-cbind(BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)

SNP<-as.matrix(row.names(data))
colnames(SNP) <- c("SNP_ukb")

data<-cbind(SNP, BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar)
rownames(data) <- c()

########################
#Save data
########################

write.csv(data, "/filepath/Percentage_GS_data_MHQ_dep_adjusted_summary_data_ukb.csv", row.names = F)


