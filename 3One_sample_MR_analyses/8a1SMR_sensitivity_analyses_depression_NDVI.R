######################################
##Script to conduct onne sample MR sensitivity analyses (IVW, MR-Egger and LADreg) for depression and NDVI
#The syntax was created by Zoe E Reed.
#The syntax was checked by Tim T. Morris
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################
library(foreign)
library(simex)
library(L1pack)
library(mr.raps)
library(coda)
library(xtable)

######################################
#Read in data
######################################

#Depression diagnosis
dep_diag<-read.csv("/filepath/all_data_dep_diag_adjusted_summary_data_ukb.csv", header=T)

#MHQ depression
dep_MHQ<-read.csv("/filepath/all_data_MHQ_dep_adjusted_summary_data_ukb.csv", header=T)

######################################
#Set up environment
######################################

set.seed(88)

#Import "ExactQ.R" function file
source("/filepath/ExactQ.R")  

#Function for generating the IGx2 
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

sign0 <- function(x){
	x[x==0] <- 1
	return(sign(x))
}

######################################
#Format data
######################################
#Depression diagnosis

#Need format: SNP, BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar

#Flip the sign for the SNP-exposure where the association is negative and then where this is the case flip the sign for the SNP-outcome association
to_flip <- sign0(dep_diag$BetaXG) == -1
dep_diag$BetaYG = dep_diag$BetaYG*sign0(dep_diag$BetaXG)
dep_diag$BetaXG = abs(dep_diag$BetaXG)

BetaXG<-dep_diag$BetaXG 
seBetaXG<-dep_diag$seBetaXG
BetaYG<-dep_diag$BetaYG 
seBetaYG<-dep_diag$seBetaYG
alphahatstar<-dep_diag$alphahatstar 
se.alphahatstar<-dep_diag$se.alphahatstar  
betastar<-dep_diag$betastar 
se.betastar<-dep_diag$se.betastar

F = BetaXG^2/seBetaXG^2

######################################
#Conduct sensitivity analyses (binary exposure)
######################################

#Implement IVW using Modified weights 
Results = weightedIVW(BetaXG,alphahatstar,seBetaXG,se.alphahatstar,tol=0.00001)
IVW_mw  = Results$RESULTS[5,1]
Qexact  = Results$QStats[4,1]
Qp      = Results$QStats[4,2]
betaIVW = betastar[1] + IVW_mw
names(Results)

#MR-Egger   
IsqGX         = Isq(BetaXG,seBetaXG)                                                   
betahat       = summary(lm(alphahatstar~ BetaXG,weights=1/se.alphahatstar^2))$coef    
MREgger       = betastar[1] + betahat[2,1]  

#Collider correction + SiMEX 
Fit           = lm(alphahatstar~BetaXG,x=TRUE,y=TRUE,weights=1/se.alphahatstar^2)        
mod.sim2      = simex(Fit,B=500,measurement.error=seBetaXG,                              
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE")
bSIMEX        = summary(mod.sim2)$coef$jackknife   
MREggersimex  = betastar[1] + bSIMEX[2,1]
EggerSE       = bSIMEX[2,2]    
Egger_intercept = bSIMEX[1,4]

#LAD regression
betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]   
LAD           = betastar[1] + betahat  

#Collider correction + SiMEX 
Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)         
mod.sim3      = simex(Fit,B=500,measurement.error=seBetaXG,       
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE",jackknife.estimation = FALSE)
bSIMEX        = mod.sim3$coef               
LADsimex      = betastar[1] + bSIMEX

#Obtain bootstrap SE for LAD regression 
Ests = NULL
for(i in 1:100){
  L     = length(BetaXG)  
  d     = sample(L,L,replace=TRUE)  
  data3 = dep_diag[d,]   
  
  BetaXG          = data3$BetaXG
  seBetaXG        = data3$seBetaXG
  alphahatstar    = data3$alphahatstar
  se.alphahatstar = data3$se.alphahatstar
  
  #LAD regression
  betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]      
  LAD           = betastar[1] + betahat                              
  
  Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)        
  mod.sim2      = simex(Fit,B=200,measurement.error=seBetaXG,
                        SIMEXvariable="BetaXG",fitting.method="quad",
                        asymptotic="FALSE",jackknife.estimation = FALSE)
  Ests[i]       = mod.sim2$coef
  print(i)               
}

seLAD = sd(Ests) 

#Sort the estimates
BetaXG          = dep_diag$BetaXG 
seBetaXG        = dep_diag$seBetaXG
alphahatstar    = dep_diag$alphahatstar
se.alphahatstar = dep_diag$se.alphahatstar

Fbar=mean(F)
Stats               = data.frame(Fbar,IsqGX,Qexact,Qp)
Estimates           = c(betaIVW,MREggersimex,LADsimex)
ColliderCorrections = Estimates-betastar[1]
SEs                 = sqrt(c(Results$RESULTS[5,2],EggerSE,seLAD)^2 + (se.betastar[1])^2)
LCIs                = Estimates - 1.96 * SEs
UCIs                = Estimates + 1.96 * SEs
pval                = 2*(1-pnorm(abs(Estimates/SEs)))

Final = data.frame(Estimates,SEs,pval,
                   row.names=c("IVW","MR-Egger","LADreg"))

fsNigx<-c(Fbar,IsqGX,NA) #Fbar for IVW and IsqGX for Egger

dep_diag_NDVI = data.frame(Estimates, LCIs, UCIs, pval,fsNigx,
                                row.names=c("IVW","MR-Egger","LADreg"))

######################################
#Format data
######################################
##MHQ
#Need format: SNP, BetaXG, seBetaXG, BetaYG, seBetaYG, alphahatstar_data, se.alphahatstar_data, betastar, se.betastar

#Flip the sign for the SNP-exposure where the association is negative and then where this is the case flip the sign for the SNP-outcome association
to_flip <- sign0(dep_MHQ$BetaXG) == -1
dep_MHQ$BetaYG = dep_MHQ$BetaYG*sign0(dep_MHQ$BetaXG)
dep_MHQ$BetaXG = abs(dep_MHQ$BetaXG)

BetaXG<-dep_MHQ$BetaXG 
seBetaXG<-dep_MHQ$seBetaXG
BetaYG<-dep_MHQ$BetaYG 
seBetaYG<-dep_MHQ$seBetaYG
alphahatstar<-dep_MHQ$alphahatstar 
se.alphahatstar<-dep_MHQ$se.alphahatstar  
betastar<-dep_MHQ$betastar 
se.betastar<-dep_MHQ$se.betastar

F = BetaXG^2/seBetaXG^2


#Implement IVW using Modified weights 
Results = weightedIVW(BetaXG,alphahatstar,seBetaXG,se.alphahatstar,tol=0.00001)
IVW_mw  = Results$RESULTS[5,1]
Qexact  = Results$QStats[4,1]
Qp      = Results$QStats[4,2]
betaIVW = betastar[1] + IVW_mw
names(Results)

#MR-Egger   
IsqGX         = Isq(BetaXG,seBetaXG)                                                   
betahat       = summary(lm(alphahatstar~ BetaXG,weights=1/se.alphahatstar^2))$coef    
MREgger       = betastar[1] + betahat[2,1]  

#Collider correction + SiMEX 
Fit           = lm(alphahatstar~BetaXG,x=TRUE,y=TRUE,weights=1/se.alphahatstar^2)        
mod.sim2      = simex(Fit,B=500,measurement.error=seBetaXG,                              
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE")
bSIMEX        = summary(mod.sim2)$coef$jackknife   
MREggersimex  = betastar[1] + bSIMEX[2,1]
EggerSE       = bSIMEX[2,2]    
Egger_intercept = bSIMEX[1,4]

#LAD regression
betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]   
LAD           = betastar[1] + betahat  

#Collider correction + SiMEX 
Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)         
mod.sim3      = simex(Fit,B=500,measurement.error=seBetaXG,       
                      SIMEXvariable="BetaXG",fitting.method="quad",
                      asymptotic="FALSE",jackknife.estimation = FALSE)
bSIMEX        = mod.sim3$coef               
LADsimex      = betastar[1] + bSIMEX

#Obtain bootstrap SE for LAD regression 
Ests = NULL
for(i in 1:100){
  L     = length(BetaXG)  
  d     = sample(L,L,replace=TRUE)  
  data3 = dep_MHQ[d,]   
  
  BetaXG          = data3$BetaXG
  seBetaXG        = data3$seBetaXG
  alphahatstar    = data3$alphahatstar
  se.alphahatstar = data3$se.alphahatstar
  
  # LAD regression
  betahat       = summary(lad(alphahatstar~-1+BetaXG))$coef[1,1]      
  LAD           = betastar[1] + betahat                              
  
  Fit           = lad(alphahatstar~-1+BetaXG,x=TRUE,y=TRUE)        
  mod.sim2      = simex(Fit,B=200,measurement.error=seBetaXG,
                        SIMEXvariable="BetaXG",fitting.method="quad",
                        asymptotic="FALSE",jackknife.estimation = FALSE)
  Ests[i]       = mod.sim2$coef
  print(i)               
}

seLAD = sd(Ests) 

#Sort the estimates
BetaXG          = dep_MHQ$BetaXG 
seBetaXG        = dep_MHQ$seBetaXG
alphahatstar    = dep_MHQ$alphahatstar
se.alphahatstar = dep_MHQ$se.alphahatstar

Fbar=mean(F)
Stats               = data.frame(Fbar,IsqGX,Qexact,Qp)
Estimates           = c(betaIVW,MREggersimex,LADsimex)
ColliderCorrections = Estimates-betastar[1]
SEs                 = sqrt(c(Results$RESULTS[5,2],EggerSE,seLAD)^2 + (se.betastar[1])^2)
LCIs                = Estimates - 1.96 * SEs
UCIs                = Estimates + 1.96 * SEs
pval                = 2*(1-pnorm(abs(Estimates/SEs)))

Final = data.frame(Estimates,SEs,pval,
                   row.names=c("IVW","MR-Egger","LADreg"))

fsNigx<-c(Fbar,IsqGX,NA) #Fbar for IVW and IsqGX for Egger

dep_MHQ_NDVI = data.frame(Estimates, LCIs, UCIs, pval,fsNigx,
                                row.names=c("IVW","MR-Egger","LADreg"))
