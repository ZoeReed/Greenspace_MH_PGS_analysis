######################################
##Script to create wellbeing PGS for MR
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
#Read in data
######################################

#Read in GWAS summary stats
#The file had a mix of tabs and spaces as seperators so this code makes sure this correctly loaded. 
#Some rows in the file had blank cells for the final column so these have been autofilled with NA
tmp<-readLines("/filepath/SWB_excl_23andMe_UKB_ldscGC_sumstats.txt")
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
well_dat$beta.exposure <- snps$Beta
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
well_dat$clumped <- FALSE #is data already clumped?
well_dat$trait <- "wellbeing" #exposure name

#Rename SNP column
well_dat <- dplyr::rename(well_dat,"SNP"="snps.SNPID")

######################################
#Generate PRS at 0.00000005
######################################

#Set working directory
setwd("/filepath/")

#Read in genome-wide significant SNP list
well_snps<-read.table("/filepath/Okbay_3_SNPs.csv", header=T)

well_5eminus8<-well_dat[(well_dat$SNP %in% well_snps$SNP), ] #3

#mrbase_grs is the function to create PGS in plink (if plink_grs=T), clumping in 500kb windows with r2 of 0.25
well_5eminus8$trait <- "well_5e8"
mrbase_grs(output = "code", exposure_dat = well_5eminus8, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_well_5e8", plink_grs = TRUE)

#Also used 0.00001 PGS previously made (for PGS analyses) as sensitivity analysis