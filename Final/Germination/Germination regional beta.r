###################################################################################################################################
#### Germination trials
#### Updared analysis May 2019


library(vegan)
library(labdsv)
library(ggplot2)
library(gridExtra)
library(indicspecies)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mydata=read.csv("Germination trials Aug 2019.csv") 


################################################################################################
# Create first matrix which includes treatments and reps
species=unique(mydata$Species...rectified)
sites <- unique(mydata$Label)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {
tdata=mydata[grep(s,mydata$Label),]
sitesp=unique(tdata$Species...rectified)
for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

tada [is.na(tada )] <- 0

write.csv(tada, file = "Germination species by site.csv") # save data out 

################################################################################################
# Create matrix based on only reps from D1 and S1

dataD1=tada[grep("D1",rownames(tada)),]
dataS1=tada[grep("S1",rownames(tada)),]
dataD1S1=rbind(dataD1,dataS1) # remove data from other replicates as different numbers of replicates undertaken


#dataD1S1=dataD1S1[!(rownames(dataD1S1) %in% c("MQ_NWW_C1_2_S1")), ]
remove= c("Spirodela spp.", "Adiantum sp.", "Cardamine flexuosa")
dataD1S1=rbind(dataD1,dataS1) # remove data from other replicates as different numbers of replicates undertaken

species=colnames(dataD1S1)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix - this adds the results for S1 and d1 together
for (s in sites) {
tdata=dataD1S1[grep(s,rownames(dataD1S1)),] # grab relevant rows
if(is.null(nrow(tdata))){
tada[grep(s,rownames(tada)),] <-tdata
} else {
tdata2=colSums(tdata) # sum together results from S1 and D1
tada[grep(s,rownames(tada)),] <-tdata2
}}



