######################################
##Script to create depression PGS for MR
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
tmp<-readLines("/filepath/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB")
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
dep_dat$beta.exposure <- log(snps$OR)
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
#Generate PRS at 0.00000005
######################################

#Set working directory
setwd("/filepath/")

#Read in genome-wide significant SNP list
dep_snps<-read.table("/filepath/Wray_44_SNPs.csv", header=T)
dep_5eminus8<-dep_dat[(dep_dat$SNP %in% dep_snps$SNP), ] #40 SNPs

#mrbase_grs is the function to create PGS in plink (if plink_grs=T), clumping in 500kb windows with r2 of 0.25
dep_5eminus8$trait <- "dep_5e8"
mrbase_grs(output = "code", exposure_dat = dep_5eminus8, clumped = FALSE, clump_kb = 500, clump_r2 = 0.25, suffix = "_dep_5e8", plink_grs = TRUE)
