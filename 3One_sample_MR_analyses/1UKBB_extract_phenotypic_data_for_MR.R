######################################
##Script to extract variables needed from main UKBB datatset for MR
#The syntax was created by Zoe E Reed.
#The syntax was checked by ....
#R Version 3.6.2
#The code below shows all manipulations and recodes with annotations

######################################
#Load libraries (use install.packages() to install)
######################################
library(data.table)

######################################
#Read in UKBB dataset
######################################

file<-fread("/filepath/data.11393.csv", data.table=F, na.strings="")

######################################
#Subset for variables of interest
######################################

non_cancer_code_varlist<-list()
icd_primary_varlist<-list()
icd_secondary_varlist<-list()
scz_MH_online_varlist<-list()

for (i in 0:28) {
	non_cancer_code_varlist[i+1]<-paste0("20002-0.", i, sep="")
}


for (i in 0:379) {
	icd_primary_varlist[i+1]<-paste0("41202-0.", i, sep="")
}

for (i in 0:434) {
	icd_secondary_varlist[i+1]<-paste0("41204-0.", i, sep="")
}

for (i in 1:16) {
	scz_MH_online_varlist[i+1]<-paste0("20544-0.", i, sep="")
}

non_cancer_code_varlist<-unlist(non_cancer_code_varlist)
icd_primary_varlist<-unlist(icd_primary_varlist)
icd_secondary_varlist<-unlist(icd_secondary_varlist)
scz_MH_online_varlist<-unlist(scz_MH_online_varlist)

vars<-c("eid","22702-0.0", "22704-0.0", "54-0.0", "21003-0.0", "22001-0.0", "20458-0.0", "20126-0.0", non_cancer_code_varlist, icd_primary_varlist, icd_secondary_varlist, scz_MH_online_varlist)
data<-file[,vars]

non_cancer_code_colnames<-list()
icd_primary_colnames<-list()
icd_secondary_colnames<-list()
scz_MH_online_colnames<-list()

for (i in 0:28) {
	non_cancer_code_colnames[i+1]<-paste0("n_20002_0_", i, sep="")
}

for (i in 0:379) {
	icd_primary_colnames[i+1]<-paste0("n_41202_0_", i, sep="")
}

for (i in 0:434) {
	icd_secondary_colnames[i+1]<-paste0("n_41204_0_", i, sep="")
}

for (i in 1:16) {
	scz_MH_online_colnames[i+1]<-paste0("n_20544_0_", i, sep="")
}

non_cancer_code_colnames<-unlist(non_cancer_code_colnames)
icd_primary_colnames<-unlist(icd_primary_colnames)
icd_secondary_colnames<-unlist(icd_secondary_colnames)
scz_MH_online_colnames<-unlist(scz_MH_online_colnames)

######################################
#Format data
######################################

colnames(data)

colnames(data)<-c("eid", "east", "north", "centre", "age", "sex", "wellbeing", "MHQ_dep", non_cancer_code_colnames, icd_primary_colnames, icd_secondary_colnames, scz_MH_online_colnames)

data$east<-as.numeric(data$east)
data$north<-as.numeric(data$north)
data$age<-as.numeric(data$age)

#For the wellbeing variable, 1 is extremely happy and 6 is extremely unhappy
#Need to remove <0
#As a higher number = more negative, we reversed the coding by subtracting from 7
data$wellbeing[data$wellbeing <0]=NA
table(data$wellbeing)
data$wellbeing<-7-data$wellbeing
table(data$wellbeing)
data$wellbeing<-ordered(data$wellbeing)
table(data$wellbeing)

#Depression MHQ
table(data$MHQ_dep)
data$MHQ_dep[data$MHQ_dep ==1]=NA
data$MHQ_dep[data$MHQ_dep ==2]=NA

data$MHQ_dep[data$MHQ_dep ==3]=1
data$MHQ_dep[data$MHQ_dep ==4]=1
data$MHQ_dep[data$MHQ_dep ==5]=1
table(data$MHQ_dep)

#Depression (non cancer codes), make any with yes for 1286 a 1 and everyone else 0. As we do not know about those with missing data we have included them as 0 too
#Do the same for schizophrenia (1289)

for (i in 0:28) {
	data[,paste0("n_20002_0_", i, sep="")][is.na(data[,paste0("n_20002_0_", i, sep="")])]=0
}

data$depression_non_cancer<-0
data$schizophrenia_non_cancer<-0

#Loop over and replace 0 with 1 where other column in loop is 1286
for (i in 0:28) {
	data$depression_non_cancer[which(data[,paste0("n_20002_0_", i, sep="")]=="1286")]<-1
}

for (i in 0:28) {
	data$schizophrenia_non_cancer[which(data[,paste0("n_20002_0_", i, sep="")]=="1289")]<-1
}

data$depression_icd_primary<-0
data$schizophrenia_icd_primary<-0

#Loop over and replace 0 with 1 where other column in loop is 1286
for (i in 0:379) {
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F32")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F320")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F321")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F322")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F323")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F328")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F329")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F33")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F330")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F331")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F332")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F333")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F334")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F338")] <-1
	data$depression_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F339")] <-1
}

table(data$depression_icd_primary)

for (i in 0:379) {
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F20")] <-1
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F200")] <-1
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F201")] <-1
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F202")] <-1
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F203")] <-1
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F204")] <-1
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F205")] <-1
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F206")] <-1
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F208")] <-1
	data$schizophrenia_icd_primary[which(data[,paste0("n_41202_0_", i, sep="")]=="F209")] <-1
}

table(data$schizophrenia_icd_primary)

#Depression (ICD secondary codes), make any with yes for F32 or F33 a 1 and everyone else 0. As we do not know about those with missing data we have included them as 0 too and for schizophrenia

data$depression_icd_secondary<-0
data$schizophrenia_icd_secondary<-0

#Loop over and replace 0 with 1 where other column in loop is 1286
for (i in 0:434) {
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F32")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F320")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F321")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F322")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F323")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F328")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F329")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F33")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F330")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F331")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F332")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F333")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F334")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F338")] <-1
	data$depression_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F339")] <-1
}

table(data$depression_icd_secondary)

for (i in 0:434) {
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F20")] <-1
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F200")] <-1
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F201")] <-1
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F202")] <-1
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F203")] <-1
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F204")] <-1
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F205")] <-1
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F206")] <-1
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F208")] <-1
	data$schizophrenia_icd_secondary[which(data[,paste0("n_41204_0_", i, sep="")]=="F209")] <-1
}

table(data$schizophrenia_icd_secondary)

#Schizophrenia/depression mental health follow-up online, make any with yes for 2 a 1 and everyone else 0. As we do not know about those with missing data we have included them as 0 too

for (i in 1:16) {
	data[,paste0("n_20544_0_", i, sep="")][is.na(data[,paste0("n_20544_0_", i, sep="")])]=0
}

data$Scz_MH_online<-0
data$Dep_MH_online<-0

#Loop over and replace 0 with 1 where other column in loop is 2
for (i in 1:16) {
	data$Scz_MH_online[which(data[,paste0("n_20544_0_", i, sep="")]=="2")]<-1
}

for (i in 1:16) {
	data$Dep_MH_online[which(data[,paste0("n_20544_0_", i, sep="")]=="11")]<-1
}

data$depression_diagnosis<-ifelse(data$depression_non_cancer == 1 | data$depression_icd_primary == 1 | data$depression_icd_secondary ==1 | data$Dep_MH_online == 1, "1", "0")

data$schizophrenia_diagnosis<-ifelse(data$schizophrenia_non_cancer == 1 | data$schizophrenia_icd_primary == 1 | data$schizophrenia_icd_secondary ==1| data$Scz_MH_online == 1, "1", "0")

#Check data
table(data$depression_non_cancer)
table(data$depression_icd_primary)
table(data$depression_icd_secondary)
table(data$Dep_MH_online)
table(data$depression_diagnosis)

table(data$schizophrenia_non_cancer)
table(data$schizophrenia_icd_primary)
table(data$schizophrenia_icd_secondary)
table(data$Scz_MH_online)
table(data$schizophrenia_diagnosis)

#Add linker file and match
linker<-read.csv("/filepath/linker_ukb21829.csv", header=T, sep=",")
colnames(linker)<-c("IID", "eid")
data2<-merge(data, linker, by="eid", all.x=T)

dim(data2) #Check dimensions are the same
head(data2)
sum(is.na(data2$IID))

#Match NDVI and percentage greenspace datasets in
load("/filepath/greenspace_biobank_final_census_20200805.RData")

greenspace<-green_data
rm(green_data)

load("/filepath/green_bump_biobank_final_census_20200805.RData")

NDVI<-green_data
rm(green_data)

#Subset green data to variables we want
greenspace<-greenspace[,c("IID", "greenspace_percent_300m")]
NDVI<-NDVI[,c("IID", "NDVI_500m_mean")]

green_data<-merge(greenspace, NDVI, by="IID", all=T)

phen<-merge(data2, green_data, by="IID", all=T)

dim(phen)
head(phen)

#Make sure withdrawn are removed
withdrawn<-read.csv("//filepath/w21829_20200204.csv", header=F)
phen<-phen[!(phen$eid %in% withdrawn$V1), ]

#Check numbers
sum(phen$schizophrenia_icd_primary==1 &phen$schizophrenia_non_cancer==1)
sum(phen$schizophrenia_icd_primary==1 &phen$schizophrenia_non_cancer==0)

sum(phen$schizophrenia_icd_secondary==1 &phen$schizophrenia_non_cancer==1)
sum(phen$schizophrenia_icd_secondary==1 &phen$schizophrenia_non_cancer==0)

sum(phen$schizophrenia_icd_primary==0 &phen$schizophrenia_non_cancer==1)
sum(phen$schizophrenia_icd_primary==0 &phen$schizophrenia_non_cancer==0)

sum(phen$schizophrenia_icd_secondary==0 &phen$schizophrenia_non_cancer==1)
sum(phen$schizophrenia_icd_secondary==0 &phen$schizophrenia_non_cancer==0)

sum(phen$schizophrenia_icd_primary==1 &phen$Scz_MH_online==1)
sum(phen$schizophrenia_icd_primary==1 &phen$Scz_MH_online==0)

sum(phen$schizophrenia_icd_secondary==1 &phen$Scz_MH_online==1)
sum(phen$schizophrenia_icd_secondary==1 &phen$Scz_MH_online==0)

sum(phen$schizophrenia_icd_primary==0 &phen$Scz_MH_online==1)
sum(phen$schizophrenia_icd_primary==0 &phen$Scz_MH_online==0)

sum(phen$schizophrenia_icd_secondary==0 &phen$Scz_MH_online==1)
sum(phen$schizophrenia_icd_secondary==0 &phen$Scz_MH_online==0)

sum(phen$schizophrenia_non_cancer==1 &phen$Scz_MH_online==1)
sum(phen$schizophrenia_non_cancer==1 &phen$Scz_MH_online==0)

sum(phen$schizophrenia_non_cancer==0 &phen$Scz_MH_online==1)
sum(phen$schizophrenia_non_cancer==0 &phen$Scz_MH_online==0)

#Subset to variables to keep
phen<-phen[,c("IID", "eid", "centre", "age", "sex", "wellbeing", "MHQ_dep", "depression_diagnosis", "schizophrenia_diagnosis", "Scz_MH_online", "schizophrenia_non_cancer", "schizophrenia_icd_primary", "schizophrenia_icd_secondary", "Dep_MH_online", "depression_non_cancer", "depression_icd_primary", "depression_icd_secondary", "greenspace_percent_300m", "NDVI_500m_mean")]

######################################
#Save data for use in next script
######################################

save(phen, file="/filepath/phenotypic_one_sample_MR.RData")
