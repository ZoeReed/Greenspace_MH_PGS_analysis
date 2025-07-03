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
tmp<-readLines("/path/GWAS/SWB_excl_23andMe_UKB_ldscGC_sumstats.txt")
tmp<-gsub(" ", "\t", tmp)
snps<-fread(text=tmp, sep="\t", header=T, fill=T)
snps<-data.frame(snps)

######################################
#Format data frame
######################################

#Create data frame with columns labelled as specified
well_dat<-data.frame(snps$SNPID) #SNP/rsid
well_dat$id.exposure <- "wellbeing" #name of exposure
well_dat$effect_allele.exposure <- snps$Effect_allele #effect allele/A1
well_dat$other_allele.exposure <- snps$Other_allele #alternative allele/A2
well_dat$eaf.exposure <- snps$EAF #EAF
well_dat$beta.exposure <- snps$Beta #effect size (beta or log(OR))
well_dat$samplesize.exposure <- NA #N
well_dat$ncase.exposure <- NA #N cases
well_dat$ncontrol.exposure <- NA #N controls
well_dat$pval.exposure <- snps$P.value #p value
well_dat$se.exposure <- snps$SE #SE
well_dat$units.exposure <- NA #units 
well_dat$exposure <- "wellbeing" #name of exposure
well_dat$mr_keep.exposure <- TRUE
well_dat$pval_origin.exposure <- "Okbay 2016" #data source
well_dat$data_source.exposure <- "Okbay 2016"
well_dat$clumped <- FALSE #indicates if data is already clumped
well_dat$trait <- "wellbeing" #exposure name

#Rename SNP column
well_dat <- dplyr::rename(well_dat, "SNP"="snps.SNPID")

######################################
#Generate PGS at different thresholds
#0.5 0.4 0.3 0.2 0.1 0.05 0.01 0.001 0.0001 0.00001 0.000001 0.00000005
######################################

#Set working directory
setwd("/path/")

#mrbase_grs is the function to create PGS in plink (if plink_grs=T), clumping in 500kb windows with r2 of 0.25
well_5<-well_dat[well_dat$pval.exposure<=0.5,]
well_5$trait <- "well_5"
mrbase_grs(output = "code",exposure_dat = well_5, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_5", plink_grs = TRUE)

well_4<-well_dat[well_dat$pval.exposure<=0.4,]
well_4$trait <- "well_4"
mrbase_grs(output = "code",exposure_dat = well_4, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_4", plink_grs = TRUE)

well_3<-well_dat[well_dat$pval.exposure<=0.3,]
well_3$trait <- "well_3"
mrbase_grs(output = "code",exposure_dat = well_3, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_3", plink_grs = TRUE)

well_2<-well_dat[well_dat$pval.exposure<=0.2,]
well_2$trait <- "well_2"
mrbase_grs(output = "code",exposure_dat = well_2, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_2", plink_grs = TRUE)

well_1<-well_dat[well_dat$pval.exposure<=0.1,]
well_1$trait <- "well_1"
mrbase_grs(output = "code",exposure_dat = well_1, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_1", plink_grs = TRUE)

well_05<-well_dat[well_dat$pval.exposure<=0.05,]
well_05$trait <- "well_05"
mrbase_grs(output = "code",exposure_dat = well_05, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_05", plink_grs = TRUE)

well_01<-well_dat[well_dat$pval.exposure<=0.01,]
well_01$trait <- "well_01"
mrbase_grs(output = "code",exposure_dat = well_01, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_01", plink_grs = TRUE)

well_001<-well_dat[well_dat$pval.exposure<=0.001,]
well_001$trait <- "well_001"
mrbase_grs(output = "code",exposure_dat = well_001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_001", plink_grs = TRUE)

well_0001<-well_dat[well_dat$pval.exposure<=0.0001,]
well_0001$trait <- "well_0001"
mrbase_grs(output = "code",exposure_dat = well_0001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_0001", plink_grs = TRUE)

well_00001<-well_dat[well_dat$pval.exposure<=0.00001,]
well_00001$trait <- "well_00001"
mrbase_grs(output = "code",exposure_dat = well_00001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_00001", plink_grs = TRUE)

well_000001<-well_dat[well_dat$pval.exposure<=0.000001,]
well_000001$trait <- "well_000001"
mrbase_grs(output = "code",exposure_dat = well_000001, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_000001", plink_grs = TRUE)

well_5eminus8<-well_dat[well_dat$pval.exposure<=5e-8,]
well_5eminus8$trait <- "well_5e8"
mrbase_grs(output = "code",exposure_dat = well_5eminus8, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_5e8", plink_grs = TRUE)
