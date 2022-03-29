######################################
##Script to create PGS for depression at different thresholds
#The syntax was created by Zoe E Reed.
#The syntax was checked by ....
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################
library("plyr")
library("dplyr")
library("data.table")
source("/filepath/mrbase_grs_v2.06.R")

######################################
#Load data
######################################

#Set working directory
setwd("/filepath")

#Read in GWAS summary stats
#The file had a mix of tabs and spaces as seperators so this code makes sure this correctly loaded. 
#Some rows in the file had blank cells for the final column so these have been autofilled with NA
tmp<-readLines("daner_pgc_mdd_meta_w2_no23andMe_rmUKBB")
tmp<-gsub(" ", "\t", tmp)
snps<-fread(text=tmp, sep="\t", header=T, fill=T)

######################################
#Format data frame
######################################

#Create data frame with columns labelled as specified
dep_dat<-data.frame(snps$SNP) #SNP/rsid
dep_dat$id.exposure <- "dep" #name of exposure
dep_dat$effect_allele.exposure <- snps$A1 #effect allele/A1
dep_dat$other_allele.exposure <- snps$A2 #alternative allele/A2
dep_dat$eaf.exposure <- snps$FRQ_A_45396 #EAF
dep_dat$beta.exposure <- log(snps$OR) #effect size (beta or log(OR))
dep_dat$samplesize.exposure <- NA #N
dep_dat$ncase.exposure <- snps$Nca #N cases
dep_dat$ncontrol.exposure <- snps$Nco #N controls
dep_dat$pval.exposure <- snps$P #p value
dep_dat$se.exposure <- snps$SE #SE
dep_dat$units.exposure <- NA #units 
dep_dat$exposure <- "dep" #name of exposure
dep_dat$mr_keep.exposure <- TRUE
dep_dat$pval_origin.exposure <- "Wray 2018" #data source
dep_dat$data_source.exposure <- "Wray 2018"
dep_dat$clumped <- FALSE #indicates if data is already clumped
dep_dat$trait <- "dep" #exposure name

#Rename SNP column
dep_dat <- dplyr::rename(dep_dat,"SNP"="snps.SNP")


######################################
#Generate PGS at different thresholds
#0.5 0.4 0.3 0.2 0.1 0.05 0.01 0.001 0.0001 0.00001 0.000001 0.00000005
######################################

#Set working directory
setwd("/filepath")

#mrbase_grs is the function to create PGS in plink (if plink_grs=T), clumping in 500kb windows with r2 of 0.25
dep_5<-dep_dat[dep_dat$pval.exposure<=0.5,]
dep_5$trait <- "dep_5"
mrbase_grs(output = "code",exposure_dat = dep_5, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_5", plink_grs = TRUE)

dep_4<-dep_dat[dep_dat$pval.exposure<=0.4,]
dep_4$trait <- "dep_4"
mrbase_grs(output = "code",exposure_dat = dep_4, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_4", plink_grs = TRUE)

dep_3<-dep_dat[dep_dat$pval.exposure<=0.3,]
dep_3$trait <- "dep_3"
mrbase_grs(output = "code",exposure_dat = dep_3, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_3", plink_grs = TRUE)

dep_2<-dep_dat[dep_dat$pval.exposure<=0.2,]
dep_2$trait <- "dep_2"
mrbase_grs(output = "code",exposure_dat = dep_2, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_2", plink_grs = TRUE)

dep_1<-dep_dat[dep_dat$pval.exposure<=0.1,]
dep_1$trait <- "dep_1"
mrbase_grs(output = "code",exposure_dat = dep_1, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_1", plink_grs = TRUE)

dep_05<-dep_dat[dep_dat$pval.exposure<=0.05,]
dep_05$trait <- "dep_05"
mrbase_grs(output = "code",exposure_dat = dep_05, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_05", plink_grs = TRUE)

dep_01<-dep_dat[dep_dat$pval.exposure<=0.01,]
dep_01$trait <- "dep_01"
mrbase_grs(output = "code",exposure_dat = dep_01, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_01", plink_grs = TRUE)

dep_001<-dep_dat[dep_dat$pval.exposure<=0.001,]
dep_001$trait <- "dep_001"
mrbase_grs(output = "code",exposure_dat = dep_001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_001", plink_grs = TRUE)

dep_0001<-dep_dat[dep_dat$pval.exposure<=0.0001,]
dep_0001$trait <- "dep_0001"
mrbase_grs(output = "code",exposure_dat = dep_0001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_0001", plink_grs = TRUE)

dep_00001<-dep_dat[dep_dat$pval.exposure<=0.00001,]
dep_00001$trait <- "dep_00001"
mrbase_grs(output = "code",exposure_dat = dep_00001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_00001", plink_grs = TRUE)

dep_000001<-dep_dat[dep_dat$pval.exposure<=0.000001,]
dep_000001$trait <- "dep_000001"
mrbase_grs(output = "code",exposure_dat = dep_000001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_000001", plink_grs = TRUE)

dep_5eminus8<-dep_dat[dep_dat$pval.exposure<=5e-8,]
dep_5eminus8$trait <- "dep_5e8"
mrbase_grs(output = "code",exposure_dat = dep_5eminus8, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_5e8", plink_grs = TRUE)
