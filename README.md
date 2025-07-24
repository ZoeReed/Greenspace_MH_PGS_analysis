# Greenspace_MH_PGS_analysis
In this study we investigated the extent to which confounding via geography may have an impact on analyses including polygenic indicators of mental health by considering within-area and contextual effects. We conducted observational and polygenic risk score analyses using UK Biobank data.

Details of how to access these data can be found in the corresponding paper.

The repository includes the following files in the following folders:

For 1Data_preparation:
Script 1 is selecting variables of interest from the main UK Biobank dataset
Script 2 is merging in greenspace measures into a single data frame (NDVI and percentage greenspace)
Script 3 is matching to the area level variables
Script 4 produces a map for a supplementary figure (map of UK Biobank participants per area)
Script 5 obtains the R squared for Principal components and produces a plot of this for a supplementary figure
Script 6 extracts the relevant phenotypic data related to depression, wellbeing and schizophrenia from the main UKBB dataset

For 2PGS_analyses:
Scripts 1a-1c create PGS for depression, wellbeing and schizophrenia
Scripts 2a-2c run the phenotypic models for the depression, wellbeing and schizophrenia phenotypes and greenspace outcomes (NDVI and percentage greenspace)
Scripts 3a-3c run the models (including multilevel models) for the depression, wellbeing and schizophrenia PGS and greenspace outcomes (NDVI and percentage greenspace). Script 3cGreenspace_and_schizophrenia_PRS_v2.R also includes code to generate a figure for the manuscript (ranked MSOA level residuals, ranked MSOA residuals from Mundlak model filtered to solely display the top and bottom 10 MSOAs and mapped locations of the extreme MSOA residuals across the dataset).

For 3Fig_generation:
Script fig_generation.R is used to create additional plots included as figures in the manuscript

PRS were created using this pipeline: https://github.com/sean-harrison-bristol/UK_Biobank_PRS/tree/master.

Syntax files can be opened as text files, but are only usable in R.

Currently (June 2025) R software can be downloaded for free here:

https://www.r-project.org/
