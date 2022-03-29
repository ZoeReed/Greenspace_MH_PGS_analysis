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
#Script to extract SNPs for risk scores
#The syntax was created by Zoe E Reed.
#The syntax was checked by Tim T. Morris
#----------------------------------

#----------------------------------
#Load required modules and set up environment
#----------------------------------

module load libs/cuda/9.0-gcc-5.4.0-2.26
module load  libs/gsl/2.5-gcc-5.5.0
module load  apps/qctool/2.0rc4

cd "/filepath"

export DIR1="/filepath/"
export DIR2="/filepath/"

#----------------------------------
#Extract data based on SNP lists
#----------------------------------

#Depression
for i in {01..22}
do
qctool -g $DIR1/data.chr"$i".bgen -og $DIR2/"$i"_Dep_all_snp.dosage.gen -s $DIR1/data.chr1-22.sample -omit-chromosome -incl-rsids $DIR2/Dep_snplist.txt -excl-samples $DIR2/Exclusions/biobank_excl_list.txt
done

rm $DIR2/Dep-snps-out.gen
touch $DIR2/Dep-snps-out.gen

for i in {01..22}
do
cat $DIR2/${i}_Dep_all_snp.dosage.gen >> $DIR2/Dep-snps-out.gen
done

#Wellbeing
for i in {01..22}
do
qctool -g $DIR1/data.chr"$i".bgen -og $DIR2/"$i"_Wellbeing_all_snp.dosage.gen -s $DIR1/data.chr1-22.sample -omit-chromosome -incl-rsids $DIR2/Wellbeing_snplist.txt -excl-samples $DIR2/Exclusions/biobank_excl_list.txt
done

rm $DIR2/Wellbeing-snps-out.gen
touch $DIR2/Wellbeing-snps-out.gen

for i in {01..22}
do
cat $DIR2/${i}_Wellbeing_all_snp.dosage.gen >> $DIR2/Wellbeing-snps-out.gen
done

#Wellbeing 1x10-05
for i in {01..22}
do
qctool -g $DIR1/data.chr"$i".bgen -og $DIR2/"$i"_Wellbeing_1x10min05_all_snp.dosage.gen -s $DIR1/data.chr1-22.sample -omit-chromosome -incl-rsids $DIR2/Wellbeing_1x10min05_snplist.txt -excl-samples $DIR2/Exclusions/biobank_excl_list.txt
done

rm $DIR2/Wellbeing-1x10min05-snps-out.gen
touch $DIR2/Wellbeing-1x10min05-snps-out.gen

for i in {01..22}
do
cat $DIR2/${i}_Wellbeing_1x10min05_all_snp.dosage.gen >> $DIR2/Wellbeing-1x10min05-snps-out.gen
done

#Schizophrenia
for i in {01..22}
do
qctool -g $DIR1/data.chr"$i".bgen -og $DIR2/"$i"_Scz_all_snp.dosage.gen -s $DIR1/data.chr1-22.sample -omit-chromosome -incl-rsids $DIR2/Scz_snplist.txt -excl-samples $DIR2/Exclusions/biobank_excl_list.txt
done

rm $DIR2/Scz-snps-out.gen
touch $DIR2/Scz-snps-out.gen

for i in {01..22}
do
cat $DIR2/${i}_Scz_all_snp.dosage.gen >> $DIR2/Scz-snps-out.gen
done

