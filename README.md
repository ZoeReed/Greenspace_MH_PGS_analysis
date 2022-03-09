# Greenspace_MH_PGS_analysis
Code to run greenspace and mental health PGS analyses

For 1Data_preparation:  
Script 1 is selecting variables of interest from the main UK Biobank dataset  
Scripts 2a and 2b are for merging in the two greenspace measures (NDVI and percentage greenspace)  
Scripts 3a and 3b are for matching to the area level variables (neighbourhood and district)  
Script 4 produces maps for supplementary figure (maps of UK Biobank participants per area)  
Script 5 obtains the R squared for Principal components and produces a plot of this for supplementary figure  

For 2PGS_analyses:  
Scripts 1a-1b create PGS for depression, wellbeing and schizophrenia  
Scripts 2a-2c run the models (including multilevel models) for the depression, wellbeing and schizophrenia PGS and greenspace outcomes (NDVI and percentage greenspace  

For 3One_sample_MR_analyses:  
Script 1 extracts the relevant phenotypic data related to depression, wellbeing and schizophrenia from the main UKBB dataset  
Scripts 2a-2c create the PGS to use in one-sample MR  
Scripts 3a-3c check the R squared for the PGS for relevant phenotypes and also includes phenotypic analyses  
Scripts 4a-4c conduct one-sample MR main analyses  
Script 5 extracts relevant SNP information for UKBB participants for use in one-sample MR sensitivity analyses  
Script 6 recodes the genotypes for these SNPs  
Scripts 7a-7c format the data needed for one-sample MR analyses  
Scripts 8a-8c conduct one-sample MR sensitivity analyses  
Script ExactQ.R is a script containing functions for use in scripts 8a-8c  
