######################################
##Script to create PGS for depression at different thresholds
#The syntax was created by Zoe E Reed.
#The syntax was checked by Gareth J. Griffith
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################

library("plyr")
library("dplyr")
library("data.table")
source("/path/mrbase_grs_v2.06.R")

######################################
#Load data
######################################

#Read in GWAS summary stats
#The file had a mix of tabs and spaces as seperators so this code makes sure this correctly loaded. 
#Some rows in the file had blank cells for the final column so these have been autofilled with NA
tmp<-readLines("/path/GWAS/daner_PGC_SCZ52_0513a")
tmp<-gsub(" ", "\t", tmp)
snps<-fread(text=tmp, sep="\t", header=T, fill=T)

######################################
#Format data frame
######################################

#Create data frame with columns labelled as specified
scz_dat<-data.frame(snps$SNP) #SNP/rsid
scz_dat$id.exposure <- "schizophrenia" #name of exposure
scz_dat$effect_allele.exposure <- snps$A1 #effect allele/A1
scz_dat$other_allele.exposure <- snps$A2 #alternative allele/A2
scz_dat$eaf.exposure <- snps$FRQ_A_35476 #EAF
scz_dat$beta.exposure <- log(snps$OR) #effect size (beta or log(OR))
scz_dat$samplesize.exposure <- NA #N
scz_dat$ncase.exposure <- 36989 #N cases
scz_dat$ncontrol.exposure <- 113075 #N controls
scz_dat$pval.exposure <- snps$P #p value
scz_dat$se.exposure <- snps$SE #SE
scz_dat$units.exposure <- NA #units 
scz_dat$exposure <- "schizophrenia" #name of exposure
scz_dat$mr_keep.exposure <- TRUE
scz_dat$pval_origin.exposure <- "Ripke 2014" #data source
scz_dat$data_source.exposure <- "Ripke 2014"
scz_dat$clumped <- FALSE #indicates if data is already clumped
scz_dat$trait <- "schizophrenia" #exposure name

#Rename SNP column
scz_dat <- dplyr::rename(scz_dat,"SNP"="snps.SNP")

######################################
#Generate PGS at different thresholds
#0.5 0.4 0.3 0.2 0.1 0.05 0.01 0.001 0.0001 0.00001 0.000001 0.00000005
######################################

#Set working directory
setwd("/path/")

#mrbase_grs is the function to create PGS in plink (if plink_grs=T), clumping in 500kb windows with r2 of 0.25
#NA's were introduced for one row so !is.na removes this row
scz_5<-scz_dat[scz_dat$pval.exposure<=0.5,]
scz_5$trait <- "scz_5"
scz_5<-scz_5[!is.na(scz_5$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_5, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_5", plink_grs = TRUE)

scz_4<-scz_dat[scz_dat$pval.exposure<=0.4,]
scz_4$trait <- "scz_4"
scz_4<-scz_4[!is.na(scz_4$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_4, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_4", plink_grs = TRUE)

scz_3<-scz_dat[scz_dat$pval.exposure<=0.3,]
scz_3$trait <- "scz_3"
scz_3<-scz_3[!is.na(scz_3$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_3, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_3", plink_grs = TRUE)

scz_2<-scz_dat[scz_dat$pval.exposure<=0.2,]
scz_2$trait <- "scz_2"
scz_2<-scz_2[!is.na(scz_2$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_2, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_2", plink_grs = TRUE)

scz_1<-scz_dat[scz_dat$pval.exposure<=0.1,]
scz_1$trait <- "scz_1"
scz_1<-scz_1[!is.na(scz_1$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_1, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_1", plink_grs = TRUE)

scz_05<-scz_dat[scz_dat$pval.exposure<=0.05,]
scz_05$trait <- "scz_05"
scz_05<-scz_05[!is.na(scz_05$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_05, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_05", plink_grs = TRUE)

scz_01<-scz_dat[scz_dat$pval.exposure<=0.01,]
scz_01$trait <- "scz_01"
scz_01<-scz_01[!is.na(scz_01$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_01, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_01", plink_grs = TRUE)

scz_001<-scz_dat[scz_dat$pval.exposure<=0.001,]
scz_001$trait <- "scz_001"
scz_001<-scz_001[!is.na(scz_001$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_001", plink_grs = TRUE)

scz_0001<-scz_dat[scz_dat$pval.exposure<=0.0001,]
scz_0001$trait <- "scz_0001"
scz_0001<-scz_0001[!is.na(scz_0001$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_0001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_0001", plink_grs = TRUE)

scz_00001<-scz_dat[scz_dat$pval.exposure<=0.00001,]
scz_00001$trait <- "scz_00001"
scz_00001<-scz_00001[!is.na(scz_00001$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_00001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_00001", plink_grs = TRUE)

scz_000001<-scz_dat[scz_dat$pval.exposure<=0.000001,]
scz_000001$trait <- "scz_000001"
scz_000001<-scz_000001[!is.na(scz_000001$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_000001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_000001", plink_grs = TRUE)

scz_5eminus8<-scz_dat[scz_dat$pval.exposure<=5e-8,]
scz_5eminus8$trait <- "scz_5e8"
scz_5eminus8<-scz_5eminus8[!is.na(scz_5eminus8$clumped),]
mrbase_grs(output = "code",exposure_dat = scz_5eminus8, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_5e8", plink_grs = TRUE)
