######################################
##Script to create schizophrenia PGS for MR
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
tmp<-readLines("/filepath/daner_PGC_SCZ52_0513a")
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
#Generate PRS at 0.00000005
######################################

#Set working directory
setwd("/filepath/")

scz_5eminus8<-scz_dat[scz_dat$pval.exposure<=5e-8,] #13,476 snps

#mrbase_grs is the function to create PGS in plink (if plink_grs=T), clumping in 500kb windows with r2 of 0.25
scz_5eminus8$trait <- "scz_5e8"
mrbase_grs(output = "code", exposure_dat = scz_5eminus8, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_scz_5e8", plink_grs = TRUE)

#Load in exposure data here and make all clumped FALSE as there was an issue with this
tmp<-read.csv("exposure_dat_scz_5e8.csv", header=T)
tmp$clumped<-"FALSE"
write.csv(tmp, "exposure_dat_scz_5e8.csv", row.names=F)
