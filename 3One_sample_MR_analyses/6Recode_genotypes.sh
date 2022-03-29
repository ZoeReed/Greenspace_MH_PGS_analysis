#!/bin/bash

#SBATCH --job-name=extract_snps
#SBATCH -o snp_out
#SBATCH -e snp_error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=serial_verylong
#SBATCH --time=300:00:00

#----------------------------------
#Script to recode SNPs to genotype
#The syntax was created by Zoe E Reed.
#The syntax was checked by ....
#----------------------------------

#----------------------------------
#Load required modules and set up environment
#----------------------------------

module load apps/plink/2.00

cd "/filepath"

export DIR1="/filepath/"
export DIR2="/filepath/"

#----------------------------------
#Recode SNPs
#----------------------------------

#Depression
plink --gen $DIR2/Dep-snps-out.gen --sample $DIR1/data.chr1-22.sample --remove $DIR2/Exclusions/biobank_excl_list.txt --allow-extra-chr --recode A --out $DIR2/Dep_snps

#Wellbeing
plink --gen $DIR2/Wellbeing-snps-out.gen --sample $DIR1/data.chr1-22.sample --remove $DIR2/Exclusions/biobank_excl_list.txt --allow-extra-chr --recode A --out $DIR2/Wellbeing_snps

#Wellbeing 1x10-05
plink --gen $DIR2/Wellbeing-1x10min05-snps-out.gen --sample $DIR1/data.chr1-22.sample --remove $DIR2/Exclusions/biobank_excl_list.txt --allow-extra-chr --recode A --out $DIR2/Wellbeing_1x10min05_snps

#Schizophrenia
plink --gen $DIR2/Scz-snps-out.gen --sample $DIR1/data.chr1-22.sample --remove $DIR2/Exclusions/biobank_excl_list.txt --allow-extra-chr --recode A --out $DIR2/Scz_snps
