######################################################################################################################################
### Script used to undertake basic analysis of wetlands with PA data across all years
### Cassie James (JCU)
### Updates July 2018

library(vegan)
library(labdsv)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(betapart)
library(indicspecies)

######################################################################################################################################
# Steps
# Loads data from different wetlands (species lists are different but codes etc should all be the same
# Creates a single species list by merging all species from all datasets and selecting unique. This is the standardised column headings 
# For each wetland creates a new dataframe using the standardised column headings and numbers but uses the individual wetlands as the rows (this will allow us to row bind later on)
# Rbinds individual datasets together for All wetland analysis
# Clean up species codes to only include simple codes with uncertain id's removed
# Initial preliminary analysis of beta diversity
# nMDS of all wetlands undertaken on P/A 
# Indicator species analysis for P/A on all wetlands
# Separated data into broad functional groups of wetland versus terrestrial and native versus exotic
# Undertaken nMDS on borad functional groups and exotic status separately
# nMDS of all wetlands by year won't converge
# nMDS of hattah and LMW together by year
# nMDS of Gunbower and KP together by year
# Stuff at the bottom is the initial runs for Hattah with data coded by time since inundation and whether the year was wet or dry ...needs sorting
#
######################################################################################################################################
### Step 1 load data from different files 

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/LMW_data_csvs/"; setwd(data.dir)
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"

data.matrix.LMW=read.csv("Spp_site_matrix LMW_WL_May 2019 COMMON ONLY.csv", row.names=1) # load data - object name is 'Output' ... :)

Output.LMW =data.matrix.LMW

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.HLWL=read.csv("Spp_site_matrix Hattah_WL_May 2019 COMMON ONLY.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.HLWL =data.matrix.HLWL


data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd(data.dir)
data.matrix.Chow=read.csv("Spp_site_matrix Chowilla_WL_May 2019 COMMON ONLY.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.Chow =data.matrix.Chow

Output.Chow=Output.Chow[,-c(5)] # remove 'NA' column


data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
data.matrix.Gun=read.csv("Spp_site_matrix Gunbower_WL_May 2019 COMMON ONLY.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.Gun =data.matrix.Gun

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/KP_data_csvs/"; setwd (data.dir) 
data.matrix.KP=read.csv("Spp_site_matrix KP_WL_May 2019 COMMON ONLY.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.KP =data.matrix.KP


data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir) 
data.matrix.Barm=read.csv("Spp_site_matrix Barmah_WL_May 2019 COMMON ONLY.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.Barm =data.matrix.Barm


##########################################################################################################################################
### Step 2
### create a single matrix of data
### create a master list of species so that all datasets have the same columns

species <- unique(c(colnames(Output.LMW),colnames(Output.HLWL),colnames(Output.Chow),colnames(Output.Gun) ,colnames(Output.KP),colnames(Output.Barm))) # create single speies list for all sites

############################################# Bring in LMW data


sites <- c(rownames(Output.LMW))
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species


tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.LMW[grep(s,rownames(Output.LMW)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.LMW  = tada


############################################# Bring in Hattah data
sites <- c(rownames(Output.HLWL))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species


tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.HLWL[grep(s,rownames(Output.HLWL)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.HLWL  = tada

############################################# Bring in Chowilla data
sites <- c(rownames(Output.Chow))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.Chow[grep(s,rownames(Output.Chow)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.Chow  = tada

wetlandsALL <-rbind(tada.HLWL, tada.LMW,tada.Chow)

wetlandsALL [is.na(wetlandsALL )] <- 0

############################################## Bring in Gunbower data

sites <- c(rownames(Output.Gun))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.Gun[grep(s,rownames(Output.Gun)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.Gun  = tada

wetlandsALL <-rbind(wetlandsALL,tada.Gun)



wetlandsALL [is.na(wetlandsALL )] <- 0

############################################## Bring in KP data

sites <- c(rownames(Output.KP))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.KP[grep(s,rownames(Output.KP)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.KP  = tada

wetlandsALL <-rbind(wetlandsALL,tada.KP)



############################################## Bring in Barmah data

sites <- c(rownames(Output.Barm))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.Barm[grep(s,rownames(Output.Barm)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.Barm  = tada

tada.Barm=as.data.frame(tada)

wetlandsALL <-rbind(wetlandsALL,tada.Barm)
wetlandsALL=ceiling(wetlandsALL)
wetlandsALL[is.na(wetlandsALL)] <- 0 

# for rare species removal...
wetlandsALL<-dropspc(wetlandsALL, 10) 

wetlandsALL [is.na(wetlandsALL )] <- 0


########################################################################################################
#### Clean up species codes and use only cleaned codes that exclude species of uncertain ID

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT","BB","CR","LP","UL","MUH","BI_","UMWC","W33","SCB","MLH","WL","WW","LLIT","PD","WWW","CSW","COX","LLIM","CWH","MIH","TWI","WOO","LLGUN","GSGUN","LGGUN","IPGUN","RLGUN","BLSGUN","FBGUN","CSGUN","LRGUN","COSGUN","PR1KP","SLKP","WHKP","BWKP","CLTKP","TLKP","PSKP","PAWKP","BCKP","PRWKP","PBKP","PJWKP","BLKP","PLLKP","WLBARM","BDEADBARM","SPBARM","TLBARM","DUCKBARM","RBSBARM","TIBBARM","BGBARM","AlgaBARM","TIOBARM","LRSBARM")

Species <-as.data.frame(colnames(wetlandsALL))

species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended

my_codes_clean=mycodes$sp_code_clean

my_codes_clean=unique(my_codes_clean[my_codes_clean!=""])# removes empty cells

wetlandsALLV2 =wetlandsALL[,which(colnames(wetlandsALL) %in% my_codes_clean)] # this removes all unknown species and most sp. unless it was the only one of its genera

wetlandsALLV2PA<-decostand(wetlandsALLV2,method="pa")
write.csv(wetlandsALLV2PA , file = "Full_TLM_dataset_PA_common_only.csv") # save data out for report Appendix 

########################################################################################################
#### MDS analysis for both wetlands together with data summed across years and then turned in PA

wetlist=c("HAT_BIT","HAT_BLT", "HAT_BOT", "HAT_BRT", "HAT_CCS", "HAT_HT", "HAT_KT", "HAT_LHAT", "HAT_MOT", "HAT_NCT", "HAT_NN", "HAT_YT","LMW_BB","LMW_CR","LMW_LP","LMW_UL","LMW_MUH","LMW_BI","LMW_UMWC","LMW_W33","LMW_SCB","LMW_MLH","LMW_WL","LMW_WW",
"CHOW_GUM","CHOW_KUL","CHOW_LLIT","CHOW_PD","CHOW_WWW","CHOW_CSW","CHOW_COX","CHOW_LLIM","CHOW_CWH","CHOW_MIH","CHOW_TWI","CHOW_WOO","CHOW_MON","CHOW_CIL","CHOW_BBW",
"GUN_LL","GUN_GS","GUN_LG","GUN_IP","GUN_RL","GUN_BLS","GUN_FB","GUN_CS","GUN_LR","GUN_COS","KP_PR1","KP_SL","KP_WH","KP_BW","KP_CLT","KP_TL","KP_PS","KP_PAW","KP_BC","KP_PRW","KP_PB","KP_PJW","KP_BL","KP_PLL",
"BARM_WL","BARM_BDEAD","BARM_SP","BARM_TL","BARM_DUCK","BARM_RBS","BARM_TIB","BARM_BG","BARM_Alga","BARM_TIO","BARM_LRS")

tada=wetlandsALLV2

tada=as.data.frame(tada)

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)
tada=tada[,colSums(tada!= 0) > 0]

tadaPA<-decostand(tada,method="pa")

result<-metaMDS((tadaPA), distance="bray", autotransform=F, k=2,trymax=100)

species.scores <- as.data.frame(scores(result, "species")) 
site.scores <- as.data.frame(scores(result, "sites")) 

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$sites=rownames(tdata)
tdata$wetlandcomplex <-tdata$sites

wetlist=c("HAT_BIT","HAT_BLT", "HAT_BOT","HAT_BRT","HAT_HT","HAT_KT", "HAT_MOT","HAT_NCT", "HAT_NN", "HAT_YT")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Hattah"

wetlist=c("LMW_BB","LMW_CR","LMW_LP","LMW_UL","LMW_MUH","LMW_BI","LMW_UMWC","LMW_W33","LMW_SCB","LMW_MLH","LMW_WL","LMW_WW")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"LMW"


wetlist=c("CHOW_GUM","CHOW_KUL","CHOW_LLIT","CHOW_PD","CHOW_WWW","CHOW_CSW","CHOW_COX","CHOW_LLIM","CHOW_CWH","CHOW_MIH","CHOW_TWI","CHOW_WOO","CHOW_MON","CHOW_BBW")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Chowilla"

wetlist=c("GUN_LL","GUN_GS","GUN_LG","GUN_IP","GUN_RL","GUN_BLS","GUN_FB","GUN_CS","GUN_LR","GUN_COS")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Gunbower"

wetlist=c("KP_PR1","KP_SL","KP_WH","KP_BW","KP_CLT","KP_TL","KP_PS","KP_PAW","KP_BC","KP_PRW","KP_PB","KP_PJW","KP_BL","KP_PLL")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"KP"

wetlist=c("BARM_WL","BARM_BDEAD","BARM_SP","BARM_TL","BARM_DUCK","BARM_RBS","BARM_TIB","BARM_BG","BARM_Alga","BARM_TIO","BARM_LRS")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Barmah"

tdata$wetlandcomplex=as.factor(tdata$wetlandcomplex)

NMDS.mean.wetlands=aggregate(tdata[,1:2],list(group=tdata$wetlandcomplex),mean)	
	
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell_wetlands <- data.frame()
for(g in levels(tdata$wetlandcomplex)){
  df_ell_wetlands <- rbind(df_ell_wetlands, cbind(as.data.frame(with(tdata[tdata$wetlandcomplex==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

setwd(image.dir)

#### Plot

png(paste(image.dir,"Wetlands together NMDS May 2019 Common species only.png",sep=''),width=12.5, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

wetlandsnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = wetlandcomplex)) +
    geom_path(data=df_ell_wetlands, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.wetlands$MDS1,y=NMDS.mean.wetlands$MDS2,label=NMDS.mean.wetlands$group)
	
wetlandsnmds
dev.off()

######################################################################################################
#### Which species are most indicative of these wetland complexes
tadaPA=as.data.frame(tadaPA)

groups=rownames(tadaPA)

wetlist=c("HAT_BIT","HAT_BLT", "HAT_BOT", "HAT_BRT", "HAT_CCS", "HAT_HT", "HAT_KT", "HAT_LHAT", "HAT_MOT", "HAT_NCT", "HAT_NN", "HAT_YT")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Hattah"

wetlist=c("LMW_BB","LMW_CR","LMW_LP","LMW_UL","LMW_MUH","LMW_BI","LMW_UMWC","LMW_W33","LMW_SCB","LMW_MLH","LMW_WL","LMW_WW")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"LMW"


wetlist=c("CHOW_GUM","CHOW_KUL","CHOW_LLIT","CHOW_PD","CHOW_WWW","CHOW_CSW","CHOW_COX","CHOW_LLIM","CHOW_CWH","CHOW_MIH","CHOW_TWI","CHOW_WOO","CHOW_MON","CHOW_CIL","CHOW_BBW")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Chowilla"

wetlist=c("GUN_LL","GUN_GS","GUN_LG","GUN_IP","GUN_RL","GUN_BLS","GUN_FB","GUN_CS","GUN_LR","GUN_COS")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Gunbower"

wetlist=c("KP_PR1","KP_SL","KP_WH","KP_BW","KP_CLT","KP_TL","KP_PS","KP_PAW","KP_BC","KP_PRW","KP_PB","KP_PJW","KP_BL","KP_PLL")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"KP"

wetlist=c("BARM_WL","BARM_BDEAD","BARM_SP","BARM_TL","BARM_DUCK","BARM_RBS","BARM_TIB","BARM_BG","BARM_Alga","BARM_TIO","BARM_LRS")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Barmah"
groups=as.factor(groups)

indval = multipatt(tadaPA, groups, control = how(nperm=999))

WetlandIndicators <-(summary(indval, alpha=0.05))


#######################################################################################################
#### Ask at broadest level are differences due to turnover or species richness

library(betapart)

tada.HLWL [is.na(tada.HLWL )] <- 0

wetlist=unique(groups)

tada = matrix(NA,nrow=7,ncol=7)#define the output matrix
rownames(tada)=c("ALL",wetlist)
colnames(tada)=c("n","BTotal", "BRepl", "BRich","BTotalVar","BReplVar","BRichVar")


tdata=tadaPA

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

mine<-beta.multi(tdata,abund=F, func="soerensen")
tada[1,1] = nrow(tdata)
tada[1,2] = mine[1]
tada[1,3] = mine[2]
tada[1,4] = mine[3]
tada[1,5] = mine[4]
tada[1,6] = mine[5]
tada[1,7] = mine[6]

tdata=tadaPA[1:12,]

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

mine<-beta.multi(tdata,abund=F, func="soerensen")
tada[2,1] = nrow(tdata)
tada[2,2] = mine[1]
tada[2,3] = mine[2]
tada[2,4] = mine[3]
tada[2,5] = mine[4]
tada[2,6] = mine[5]
tada[2,7] = mine[6]

tdata=tadaPA[13:24,]

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

mine<-beta.multi(tdata,abund=F, func="soerensen")
tada[3,1] = nrow(tdata)
tada[3,2] = mine[1]
tada[3,3] = mine[2]
tada[3,4] = mine[3]
tada[3,5] = mine[4]
tada[3,6] = mine[5]
tada[3,7] = mine[6]

tdata=tadaPA[25:34,]

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

mine<-beta.multi(tdata,abund=F, func="soerensen")
tada[4,1] = nrow(tdata)
tada[4,2] = mine[1]
tada[4,3] = mine[2]
tada[4,4] = mine[3]
tada[4,5] = mine[4]
tada[4,6] = mine[5]
tada[4,7] = mine[6]

tdata=tadaPA[35:44,]

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

mine<-beta.multi(tdata,abund=F, func="soerensen")
tada[5,1] = nrow(tdata)
tada[5,2] = mine[1]
tada[5,3] = mine[2]
tada[5,4] = mine[3]
tada[5,5] = mine[4]
tada[5,6] = mine[5]
tada[5,7] = mine[6]

tdata=tadaPA[45:58,]

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

mine<-beta.multi(tdata,abund=F, func="soerensen")
tada[6,1] = nrow(tdata)
tada[6,2] = mine[1]
tada[6,3] = mine[2]
tada[6,4] = mine[3]
tada[6,5] = mine[4]
tada[6,6] = mine[5]
tada[6,7] = mine[6]

tdata=tadaPA[59:69,]

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]

mine<-beta.multi(tdata,abund=F, func="soerensen")
tada[7,1] = nrow(tdata)
tada[7,2] = mine[1]
tada[7,3] = mine[2]
tada[7,4] = mine[3]
tada[7,5] = mine[4]
tada[7,6] = mine[5]
tada[7,7] = mine[6]


write.csv(tada , file = "Wetlands beta analysis PAdata.csv") # save data out
#########################################################################################################
#### Separate data into broad functional groups and exotic status

species.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended
mycodes=mycodes[!duplicated(mycodes[,c('sp_code_clean')]),] # remove duplicates

myspeciers=as.data.frame(colnames(wetlandsALLV2))
colnames(myspeciers)="Species"
fgrps=merge(myspeciers,mycodes,by.x="Species",by.y="sp_code_simple",all.x=TRUE)

fgrps=fgrps[,c("Species","Final_PFG_allocation","Weed.status")]
fgrps$Species=as.character(fgrps$Species)

twetlandALL=as.data.frame(t(wetlandsALLV2))

fgrpswetlands=merge(fgrps,twetlandALL,by.x="Species",by.y="row.names",all.y=TRUE)

Amph=c("ATw","ATl","ARp","ARf","ATe","S","Se")
Terre=c("Tdr","Tdr")

Wetgroup=fgrpswetlands[which(fgrpswetlands$Final_PFG_allocation %in% Amph),]
Drygroup=fgrpswetlands[which(fgrpswetlands$Final_PFG_allocation %in% Terre),]
WetgroupNats=Wetgroup[which(Wetgroup$Weed.status == "FALSE"),]
WetgroupExotic=Wetgroup[which(Wetgroup$Weed.status == "TRUE"),]
DrygroupNats=Drygroup[which(Drygroup$Weed.status == "FALSE"),]
DrygroupExotic=Drygroup[which(Drygroup$Weed.status == "TRUE"),]

rownames(WetgroupNats)=WetgroupNats$Species
rownames(WetgroupExotic)=WetgroupExotic$Species
rownames(DrygroupNats)=DrygroupNats$Species
rownames(DrygroupExotic)=DrygroupExotic$Species

WetgroupNats=WetgroupNats[,-which(names(WetgroupNats) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
WetgroupExotic=WetgroupExotic[,-which(names(WetgroupExotic) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
DrygroupNats=DrygroupNats[,-which(names(DrygroupNats) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
DrygroupExotic=DrygroupExotic[,-which(names(DrygroupExotic) %in% c("Species","Final_PFG_allocation", "Weed.status"))]

WetgroupNats=as.data.frame(t(WetgroupNats))
WetgroupExotic=as.data.frame(t(WetgroupExotic))
DrygroupNats=as.data.frame(t(DrygroupNats))
DrygroupExotic=as.data.frame(t(DrygroupExotic))


#########################################################################################################
#### Separate data into broad functional groups and exotic status

tada=WetgroupNats

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)
tada=tada[,colSums(tada!= 0) > 0]

tadaPA<-decostand(tada,method="pa")
result<-metaMDS((tadaPA), distance="bray", autotransform=F, k=2,trymax=100)

species.scores <- as.data.frame(scores(result, "species")) 
site.scores <- as.data.frame(scores(result, "sites")) 

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$sites=rownames(tdata)
tdata$wetlandcomplex <-tdata$sites

wetlist=c("HAT_BIT","HAT_BLT", "HAT_BOT", "HAT_BRT", "HAT_CCS", "HAT_HT", "HAT_KT", "HAT_LHAT", "HAT_MOT", "HAT_NCT", "HAT_NN", "HAT_YT")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Hattah"

wetlist=c("LMW_BB","LMW_CR","LMW_LP","LMW_UL","LMW_MUH","LMW_BI","LMW_UMWC","LMW_W33","LMW_SCB","LMW_MLH","LMW_WL","LMW_WW")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"LMW"


wetlist=c("CHOW_GUM","CHOW_KUL","CHOW_LLIT","CHOW_PD","CHOW_WWW","CHOW_CSW","CHOW_COX","CHOW_LLIM","CHOW_CWH","CHOW_MIH","CHOW_TWI","CHOW_WOO","CHOW_MON","CHOW_CIL","CHOW_BBW")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Chowilla"

wetlist=c("GUN_LL","GUN_GS","GUN_LG","GUN_IP","GUN_RL","GUN_BLS","GUN_FB","GUN_CS","GUN_LR","GUN_COS")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Gunbower"

wetlist=c("KP_PR1","KP_SL","KP_WH","KP_BW","KP_CLT","KP_TL","KP_PS","KP_PAW","KP_BC","KP_PRW","KP_PB","KP_PJW","KP_BL","KP_PLL")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"KP"

wetlist=c("BARM_WL","BARM_BDEAD","BARM_SP","BARM_TL","BARM_DUCK","BARM_RBS","BARM_TIB","BARM_BG","BARM_Alga","BARM_TIO","BARM_LRS")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Barmah"

tdata$wetlandcomplex=as.factor(tdata$wetlandcomplex)

NMDS.mean.wetlands=aggregate(tdata[,1:2],list(group=tdata$wetlandcomplex),mean)	
	
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell_wetlands <- data.frame()
for(g in levels(tdata$wetlandcomplex)){
  df_ell_wetlands <- rbind(df_ell_wetlands, cbind(as.data.frame(with(tdata[tdata$wetlandcomplex==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}



# basic plot

png(paste(image.dir,"Wetlands Wet Natives NMDS.png",sep=''),width=12.5, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

wetlandsnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = wetlandcomplex)) +
    geom_path(data=df_ell_wetlands, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.wetlands$MDS1,y=NMDS.mean.wetlands$MDS2,label=NMDS.mean.wetlands$group)
	
wetlandsnmds
dev.off()

######################################################################################################
# Indicator species for functional group analysis

tadaPA=as.data.frame(tadaPA)

groups=rownames(tadaPA)

wetlist=c("HAT_BIT","HAT_BLT", "HAT_BOT", "HAT_BRT", "HAT_CCS", "HAT_HT", "HAT_KT", "HAT_LHAT", "HAT_MOT", "HAT_NCT", "HAT_NN", "HAT_YT")
groups[which(groups %in% wetlist)] <-"Hattah"

wetlist=c("LMW_BB","LMW_CR","LMW_LP","LMW_UL","LMW_MUH","LMW_BI","LMW_UMWC","LMW_W33","LMW_SCB","LMW_MLH","LMW_WL","LMW_WW")
groups[which(groups %in% wetlist)] <-"LMW"

wetlist=c("CHOW_GUM","CHOW_KUL","CHOW_LLIT","CHOW_PD","CHOW_WWW","CHOW_CSW","CHOW_COX","CHOW_LLIM","CHOW_CWH","CHOW_MIH","CHOW_TWI","CHOW_WOO","CHOW_MON","CHOW_CIL","CHOW_BBW")
groups[which(groups %in% wetlist)] <-"Chowilla"

wetlist=c("GUN_LL","GUN_GS","GUN_LG","GUN_IP","GUN_RL","GUN_BLS","GUN_FB","GUN_CS","GUN_LR","GUN_COS")
groups[which(groups %in% wetlist)] <-"Gunbower"

wetlist=c("KP_PR1","KP_SL","KP_WH","KP_BW","KP_CLT","KP_TL","KP_PS","KP_PAW","KP_BC","KP_PRW","KP_PB","KP_PJW","KP_BL","KP_PLL")
groups[which(groups %in% wetlist)] <-"KP"

wetlist=c("BARM_WL","BARM_BDEAD","BARM_SP","BARM_TL","BARM_DUCK","BARM_RBS","BARM_TIB","BARM_BG","BARM_Alga","BARM_TIO","BARM_LRS")
groups[which(groups %in% wetlist)] <-"Barmah"

groups=as.factor(groups)

indval = multipatt(tadaPA, groups, control = how(nperm=999))

WetlandIndicators <-(summary(indval, alpha=0.05))


#########################################################################################################
#### MDS analysis for both wetlands together with years data separated
#### Import data

species <- unique(c(colnames(Output.LMW),colnames(Output.HLWL),colnames(Output.Chow),colnames(Output.Gun) ,colnames(Output.KP))) # create single speies list for all sites
image.dir="C:/Users/jc246980/Documents/MD Vegetation/Plots/"

data.dir="C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd(data.dir)
data.matrix.LMW=read.csv("Spp_site_matrix summarised to wetland_LMW_WL_June 2018.csv", row.names=1) # load data - object name is 'Output' ... :)
Output.LMW =data.matrix.LMW

data.dir="C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.HLWL=read.csv("Spp_site_matrix summarised to wetland_HTH_WL_June 2018.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.HLWL =data.matrix.HLWL

data.dir="C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd(data.dir)
data.matrix.Chow=read.csv("Spp_site_year matrix wetland_Chowilla_WL_June 2018.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.Chow =data.matrix.Chow

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
data.matrix.Gun=read.csv("Spp_site_year matrix wetland_Gunbower_WL_June 2018.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.Gun =data.matrix.Gun

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/KP_data_csvs/"; setwd (data.dir) 
data.matrix.KP=read.csv("Spp_site_year matrix wetland_KP_WL_June 2018.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.KP=data.matrix.KP

#### Create matrix with same species names as the matrices for the different sites currently have different species lists

species <- unique(c(colnames(Output.LMW),colnames(Output.HLWL),colnames(Output.Chow),colnames(Output.Gun) ,colnames(Output.KP))) # create single speies list for all sites
sites <- c(rownames(Output.LMW))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=Output.LMW[grep(s,rownames(Output.LMW)),]

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata[1]
}}

tada.lmw  = tada

sites <- c(rownames(Output.HLWL))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.HLWL[grep(s,rownames(Output.HLWL)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.HLWL  = tada


sites <- c(rownames(Output.Chow))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.Chow[grep(s,rownames(Output.Chow)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.Chow  = tada

sites <- c(rownames(Output.Gun))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.Gun[grep(s,rownames(Output.Gun)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.Gun  = tada

sites <- c(rownames(Output.KP))

tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for(s in sites) { # 
tdata=(as.data.frame(Output.KP[grep(s,rownames(Output.KP)),]))

for (i in 1:ncol(tdata)) {
ttdata=tdata[,i]
soi=colnames(tdata)[i]
tada[grep(s,rownames(tada)),grep(soi,colnames(tada))] <-ttdata
}}

tada.KP  = tada


##########################################################################################################################
### Analysis of Hattah and LMW together - nMDS won't converge when any other data is incorporated in

wetlandsALL <-rbind(tada.HLWL,tada.lmw)

wetlandsALL[is.na(wetlandsALL)] <- 0 # replace nas with zeros in species matrix
wetlandsALL=wetlandsALL[rowSums(wetlandsALL!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)
wetlandsALL=wetlandsALL[,colSums(wetlandsALL!= 0) >0]

wetlandsALLPA<-decostand(wetlandsALL,method="pa")
#wetlandsALLPA<-dropspc(wetlandsALLPA, 2)   # remove species that occur in = < two sites/date combos
wetlandsALLPA=wetlandsALLPA[rowSums(wetlandsALLPA!= 0) > 0,]	# remove sites with no records 
wetlandsALLPA=wetlandsALLPA[,colSums(wetlandsALLPA!= 0) > 0]	# remove sites with no re
#wetlandsALLPA <-decostand(wetlandsALLPA,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative value
result<-metaMDS((wetlandsALLPA), autotransform=F, k=3,trymax=20)

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2], MDS3=result$points[,3])
tdata$group=rownames(tdata)

tdata$group[grepl("_08",tdata$group)] <- "2008"
tdata$group[grepl("_09",tdata$group)] <- "2009"
tdata$group[grepl("_10",tdata$group)] <- "2010"
tdata$group[grepl("_11",tdata$group)] <- "2011"
tdata$group[grepl("_12",tdata$group)] <- "2012"
tdata$group[grepl("_13",tdata$group)] <- "2013"
tdata$group[grepl("_14",tdata$group)] <- "2014"
tdata$group[grepl("_16",tdata$group)] <- "2016"
tdata$group[grepl("2005",tdata$group)] <- "2005"
tdata$group[grepl("2006",tdata$group)] <- "2006"
tdata$group[grepl("2007",tdata$group)] <- "2007"
tdata$group[grepl("2008",tdata$group)] <- "2008"
tdata$group[grepl("2009",tdata$group)] <- "2009"
tdata$group[grepl("2010",tdata$group)] <- "2010"
tdata$group[grepl("2011",tdata$group)] <- "2011"
tdata$group[grepl("2012",tdata$group)] <- "2012"
tdata$group[grepl("2013",tdata$group)] <- "2013"
tdata$group[grepl("2014",tdata$group)] <- "2014"
tdata$group[grepl("2015",tdata$group)] <- "2014"
tdata$group[grepl("2016",tdata$group)] <- "2016"
tdata$sites<- rownames(tdata)
tdata$sites[grepl("BB",tdata$sites)] <- "BB"
tdata$sites[grepl("CR",tdata$sites)] <- "CR"
tdata$sites[grepl("WW",tdata$sites)] <- "WW"
tdata$sites[grepl("LP",tdata$sites)] <- "LP"
tdata$sites[grepl("MUH",tdata$sites)] <- "MUH"
tdata$sites[grepl("UL",tdata$sites)] <- "UL"
tdata$sites[grepl("UMWC",tdata$sites)] <- "UMWC"
tdata$sites[grepl("MLH",tdata$sites)] <- "MLH"
tdata$sites[grepl("BI_",tdata$sites)] <- "BI"
tdata$sites[grepl("W33",tdata$sites)] <- "W33"
tdata$sites[grepl("SCB",tdata$sites)] <- "SCB"
tdata$sites[grepl("WL",tdata$sites)] <- "WL"
tdata$sites[grepl("BIT",tdata$sites)] <- "BIT"
tdata$sites[grepl("BLT",tdata$sites)] <- "BLT"
tdata$sites[grepl("BOT",tdata$sites)] <- "BOT"
tdata$sites[grepl("BRT",tdata$sites)] <- "BRT"
tdata$sites[grepl("CCS",tdata$sites)] <- "CCS"
tdata$sites[grepl("HT",tdata$sites)] <- "HT"
tdata$sites[grepl("LHAT",tdata$sites)] <- "LHAT"
tdata$sites[grepl("KT",tdata$sites)] <- "KT"
tdata$sites[grepl("YT",tdata$sites)] <- "YT"
tdata$sites[grepl("MOT",tdata$sites)] <- "MOT"
tdata$sites[grepl("NCT",tdata$sites)] <- "NCT"
tdata$sites[grepl("NN",tdata$sites)] <- "NN"

tdata$wetlandcomplex <-tdata$sites
tdata$group=as.factor(tdata$group)
tdata$sites=as.factor(tdata$sites)

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Hattah"

wetlist=c("BB","CR","LP","UL","MUH","BI","UMWC","W33","SCB","MLH","WL","WW")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"LMW"

tdata$group=as.factor(tdata$group)
tdata$sites=as.factor(tdata$sites)
tdata$wetlandcomplex=as.factor(tdata$wetlandcomplex)


NMDS.mean=aggregate(tdata[,1:2],list(group=tdata$group),mean)
NMDS.mean.sites=aggregate(tdata[,1:2],list(group=tdata$sites),mean)	
NMDS.mean.wetlands=aggregate(tdata[,1:2],list(group=tdata$wetlandcomplex),mean)	
	
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
				   
df_ell <- data.frame()
for(g in levels(tdata$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(tdata[tdata$group==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_sites <- data.frame()
for(g in levels(tdata$sites)){
  df_ell_sites <- rbind(df_ell_sites, cbind(as.data.frame(with(tdata[tdata$sites==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_wetlands <- data.frame()
for(g in levels(tdata$wetlandcomplex)){
  df_ell_wetlands <- rbind(df_ell_wetlands, cbind(as.data.frame(with(tdata[tdata$wetlandcomplex==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

png(paste(image.dir,"HTLW and LMW together NMDS wetland complex.png",sep=''),width=12.5, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

wetlandsnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = wetlandcomplex)) +
    geom_path(data=df_ell_wetlands, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.wetlands$MDS1,y=NMDS.mean.wetlands$MDS2,label=NMDS.mean.wetlands$group)
	
wetlandsnmds
dev.off()

png(paste(image.dir,"HTLW and LMW together NMDS.png",sep=''),width=25, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two column

yearnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)

sitenmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = sites)) +
    geom_path(data=df_ell_sites, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.sites$MDS1,y=NMDS.mean.sites$MDS2,label=NMDS.mean.sites$group)
	
	grid.arrange(yearnmds,sitenmds, ncol = 2)
	
dev.off()

png(paste(image.dir,"ALL Wetland NMDS sites separate.png",sep=''),width=20, height=60, units='cm', res=300, pointsize=20, bg='white')
		
ggplot(data = tdata, aes(MDS1, MDS2)) +
geom_point(data=tdata[,c("MDS1","MDS2")],size=1,aes(x=MDS1,y=MDS2),colour="grey")+
geom_point(aes(color = group),size=2)+
facet_wrap(~sites,ncol=2)+geom_path(linetype = "dashed",colour="darkgrey") 

dev.off()

##########################################################################################################################
### Analysis of Gunbower and KP together - nMDS won't converge when any other data is incorporated in

wetlandsALL <-rbind(tada.Gun,tada.KP)


wetlandsALL[is.na(wetlandsALL)] <- 0 # replace nas with zeros in species matrix
wetlandsALL=wetlandsALL[rowSums(wetlandsALL!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)
wetlandsALL=wetlandsALL[,colSums(wetlandsALL!= 0) >0]

wetlandsALLPA<-decostand(wetlandsALL,method="pa")
#wetlandsALLPA<-dropspc(wetlandsALLPA, 2)   # remove species that occur in = < two sites/date combos
wetlandsALLPA=wetlandsALLPA[rowSums(wetlandsALLPA!= 0) > 0,]	# remove sites with no records 
wetlandsALLPA=wetlandsALLPA[,colSums(wetlandsALLPA!= 0) > 0]	# remove sites with no re
#wetlandsALLPA <-decostand(wetlandsALLPA,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative value
result<-metaMDS((wetlandsALLPA), autotransform=F, k=3,trymax=20)


tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2], MDS3=result$points[,3])
tdata$group=rownames(tdata)

tdata$group[grepl("2004",tdata$group)] <- "2004"
tdata$group[grepl("2005",tdata$group)] <- "2005"
tdata$group[grepl("2006",tdata$group)] <- "2006"
tdata$group[grepl("2007",tdata$group)] <- "2007"
tdata$group[grepl("2008",tdata$group)] <- "2008"
tdata$group[grepl("2009",tdata$group)] <- "2009"
tdata$group[grepl("2010",tdata$group)] <- "2010"
tdata$group[grepl("2011",tdata$group)] <- "2011"
tdata$group[grepl("2012",tdata$group)] <- "2012"
tdata$group[grepl("2013",tdata$group)] <- "2013"
tdata$group[grepl("2014",tdata$group)] <- "2014"
tdata$group[grepl("2015",tdata$group)] <- "2015"
tdata$group[grepl("2016",tdata$group)] <- "2016"


tdata$sites<- rownames(tdata)
tdata$sites[grepl("LL_",tdata$sites)] <- "LL"
tdata$sites[grepl("GS",tdata$sites)] <- "GS"
tdata$sites[grepl("LG",tdata$sites)] <- "LG"
tdata$sites[grepl("IP",tdata$sites)] <- "IP"
tdata$sites[grepl("RL",tdata$sites)] <- "RL"
tdata$sites[grepl("BLS",tdata$sites)] <- "BLS"
tdata$sites[grepl("FB",tdata$sites)] <- "FB"
tdata$sites[grepl("CS",tdata$sites)] <- "CS"
tdata$sites[grepl("LR",tdata$sites)] <- "LR"
tdata$sites[grepl("COS",tdata$sites)] <- "COS"
tdata$sites[grepl("PR1",tdata$sites)] <- "PR1"
tdata$sites[grepl("SL",tdata$sites)] <- "SL"
tdata$sites[grepl("WH",tdata$sites)] <- "WH"
tdata$sites[grepl("BW",tdata$sites)] <- "BW"
tdata$sites[grepl("CLT",tdata$sites)] <- "CLT"
tdata$sites[grepl("TL",tdata$sites)] <- "TL"
tdata$sites[grepl("PS",tdata$sites)] <- "PS"
tdata$sites[grepl("PAW",tdata$sites)] <- "PAW"
tdata$sites[grepl("BC",tdata$sites)] <- "BC"
tdata$sites[grepl("PRW",tdata$sites)] <- "PRW"
tdata$sites[grepl("PB",tdata$sites)] <- "PB"
tdata$sites[grepl("PJW",tdata$sites)] <- "PJW"
tdata$sites[grepl("BL_",tdata$sites)] <- "BL"
tdata$sites[grepl("PLL",tdata$sites)] <- "PLL"

tdata$wetlandcomplex <-tdata$sites
tdata$group=as.factor(tdata$group)
tdata$sites=as.factor(tdata$sites)

GUNsites=c("LL","GS","LG","IP","RL","BLS","FB","CS","LR","COS")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% GUNsites)] <-"Gunbower"
KPsites=c("PR1","SL","WH","BW","CLT","TL","PS","PAW","BC","PRW","PB","PJW","BL","PLL")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% KPsites)] <-"KP"

tdata$group=as.factor(tdata$group)
tdata$sites=as.factor(tdata$sites)
tdata$wetlandcomplex=as.factor(tdata$wetlandcomplex)


NMDS.mean=aggregate(tdata[,1:2],list(group=tdata$group),mean)
NMDS.mean.sites=aggregate(tdata[,1:2],list(group=tdata$sites),mean)	
NMDS.mean.wetlands=aggregate(tdata[,1:2],list(group=tdata$wetlandcomplex),mean)	
	
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
				   
df_ell <- data.frame()
for(g in levels(tdata$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(tdata[tdata$group==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_sites <- data.frame()
for(g in levels(tdata$sites)){
  df_ell_sites <- rbind(df_ell_sites, cbind(as.data.frame(with(tdata[tdata$sites==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_wetlands <- data.frame()
for(g in levels(tdata$wetlandcomplex)){
  df_ell_wetlands <- rbind(df_ell_wetlands, cbind(as.data.frame(with(tdata[tdata$wetlandcomplex==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

png(paste(image.dir,"Gun andKP together NMDS wetland complex.png",sep=''),width=12.5, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

wetlandsnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = wetlandcomplex)) +
    geom_path(data=df_ell_wetlands, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.wetlands$MDS1,y=NMDS.mean.wetlands$MDS2,label=NMDS.mean.wetlands$group)
	
wetlandsnmds
dev.off()

png(paste(image.dir,"Gun and KP together NMDS.png",sep=''),width=25, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two column

yearnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)

sitenmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = sites)) +
    geom_path(data=df_ell_sites, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.sites$MDS1,y=NMDS.mean.sites$MDS2,label=NMDS.mean.sites$group)
	
	grid.arrange(yearnmds,sitenmds, ncol = 2)
	
dev.off()

tdataGun=tdata[which(tdata$wetlandcomplex=="Gunbower"),]

png(paste(image.dir,"Gunbower Wetland NMDS sites separate based on KP and GUN nMDS.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
		
ggplot(data = tdataGun, aes(MDS1, MDS2)) +
geom_point(data=tdata[,c("MDS1","MDS2")],size=1,aes(x=MDS1,y=MDS2),colour="grey")+
geom_point(aes(color = group),size=2)+
facet_wrap(~sites,ncol=2)+geom_path(linetype = "dashed",colour="darkgrey") 

dev.off()


tdataKP=tdata[which(tdata$wetlandcomplex=="KP"),]

png(paste(image.dir,"KP Wetland NMDS sites separate based on KP and GUN nMDS.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
		
ggplot(data = tdataKP, aes(MDS1, MDS2)) +
geom_point(data=tdata[,c("MDS1","MDS2")],size=1,aes(x=MDS1,y=MDS2),colour="grey")+
geom_point(aes(color = group),size=2)+
facet_wrap(~sites,ncol=2)+geom_path(linetype = "dashed",colour="darkgrey") 

dev.off()



########################################################################################################
#### MDS analysis for Hattah, LMW and Chowilla only
tada=wetlandsALL
wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT","BB","CR","LP","UL","MUH","BI_","UMWC","W33","SCB","MLH","WL","WW","LLIT","PD","WWW","CSW","COX","LLIM","CWH","MIH","TWI","WOO")


tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)
tada=tada[,colSums(tada!= 0) > 0]

tadaPA<-decostand(tada,method="pa")
result<-metaMDS((tadaPA), distance="bray", autotransform=F, k=2,trymax=100)

species.scores <- as.data.frame(scores(result, "species")) 
site.scores <- as.data.frame(scores(result, "sites")) 

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$sites=rownames(tdata)
tdata$wetlandcomplex <-tdata$sites

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Hattah"

wetlist=c("BB","CR","LP","UL","MUH","BI_","UMWC","W33","SCB","MLH","WL","WW")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"LMW"

wetlist=c("LLIT","PD","WWW","CSW","COX","LLIM","CWH","MIH","TWI","WOO")
tdata$wetlandcomplex[which(tdata$wetlandcomplex %in% wetlist)] <-"Chowilla"


tdata$wetlandcomplex=as.factor(tdata$wetlandcomplex)

NMDS.mean.wetlands=aggregate(tdata[,1:2],list(group=tdata$wetlandcomplex),mean)	
	
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell_wetlands <- data.frame()
for(g in levels(tdata$wetlandcomplex)){
  df_ell_wetlands <- rbind(df_ell_wetlands, cbind(as.data.frame(with(tdata[tdata$wetlandcomplex==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

png(paste(image.dir,"Wetlands Hattah LMW and Chowilla together NMDS PA all years combined.png",sep=''),width=12.5, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

wetlandsnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = wetlandcomplex)) +
    geom_path(data=df_ell_wetlands, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.wetlands$MDS1,y=NMDS.mean.wetlands$MDS2,label=NMDS.mean.wetlands$group)+
	geom_text(data=site.scores,aes(x=NMDS1,y=NMDS2,label=rownames(site.scores)),alpha=0.5)+ guides(shape=FALSE)
	
wetlandsnmds
dev.off()
 # better plot without site labels which overlap too much

png(paste(image.dir,"Wetlands together NMDS.png",sep=''),width=12.5, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

wetlandsnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = wetlandcomplex)) +
    geom_path(data=df_ell_wetlands, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.wetlands$MDS1,y=NMDS.mean.wetlands$MDS2,label=NMDS.mean.wetlands$group)
	
wetlandsnmds
dev.off()


#####################################################################################################################################################
wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT","BB","CR","LP","UL","MUH","BI_","UMWC","W33","SCB","MLH","WL","WW")

tada = matrix(NA,nrow=length(wetlist),ncol=5)#define the output matrix
rownames(tada)=wetlist
colnames(tada)=c("n","BSor", "BSim", "BNes","BTbc")

right = function(text, num_char) {
  substr(text, nchar(text) - (num_char-1), nchar(text))
}
		
for(w in wetlist) { # 
tdata=wetlandsALL[grep(w,rownames(wetlandsALL)),]

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

Baselga<-nestedbetasor(tdata)
tada[grep(w,rownames(tada)),1] = nrow(tdata)
tada[grep(w,rownames(tada)),2] = Baselga[[3]]
tada[grep(w,rownames(tada)),3] = Baselga[[1]]
tada[grep(w,rownames(tada)),4] = Baselga[[2]]
tada[grep(w,rownames(tada)),5] = mean(vegdist(tdata))
}

write.csv(tada , file = "All wetlands_betawithinsites.csv") # save data out

######################################################################################################################################


library(BAT)

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT","BB","CR","LP","UL","MUH","BI_","UMWC","W33","SCB","MLH","WL","WW")

tada = matrix(NA,nrow=length(wetlist),ncol=7)#define the output matrix
rownames(tada)=wetlist
colnames(tada)=c("n","BTotal", "BRepl", "BRich","BTotalVar","BReplVar","BRichVar")

right = function(text, num_char) {
  substr(text, nchar(text) - (num_char-1), nchar(text))
}
		
for(w in wetlist) { # 
tdata=wetlandsALL[grep(w,rownames(wetlandsALL)),]

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

mine<-beta.multi(tdata,abund=F, func="soerensen")
tada[grep(w,rownames(tada)),1] = nrow(tdata)
tada[grep(w,rownames(tada)),2] = mine[1]
tada[grep(w,rownames(tada)),3] = mine[2]
tada[grep(w,rownames(tada)),4] = mine[3]
tada[grep(w,rownames(tada)),5] = mine[4]
tada[grep(w,rownames(tada)),6] = mine[5]
tada[grep(w,rownames(tada)),7] = mine[6]
}

write.csv(tada , file = "All wetlands_betawithinsites.csv") # save data out


data.dir="C:/Users/jc246980/Documents/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd(data.dir)
Hydrodata=data.frame(read.csv("Hydraulics_HTH_WL.csv")) # load in hydraulic history for each site


Hydrodata$wetyear <-"FALSE"
Hydrodata$wetyear[which(Hydrodata$WetYear>0)]<-"TRUE"

for(w in wetlist) { # 

if(w=="BIT"){
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="BRT"){
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="BOT"){
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="BLT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
}
if(w=="CCS"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2011",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2012",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}

if(w=="NCT"){
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
}
if(w=="HT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2012",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="KT"){
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2014",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
}
if(w=="LHAT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
}
if(w=="MOT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="NN"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2011",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2012",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2014",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="YT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2011",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2012",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2014",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
}

Hydrodata$sites<- as.character(Hydrodata$Site.ID)


Hydrodata$sites[grepl("BIT",Hydrodata$sites)] <- "BIT"
Hydrodata$sites[grepl("BLT",Hydrodata$sites)] <- "BLT"
Hydrodata$sites[grepl("BOT",Hydrodata$sites)] <- "BOT"
Hydrodata$sites[grepl("BRT",Hydrodata$sites)] <- "BRT"
Hydrodata$sites[grepl("CCS",Hydrodata$sites)] <- "CCS"
Hydrodata$sites[grepl("HT",Hydrodata$sites)] <- "HT"
Hydrodata$sites[grepl("LHAT",Hydrodata$sites)] <- "LHAT"
Hydrodata$sites[grepl("KT",Hydrodata$sites)] <- "KT"
Hydrodata$sites[grepl("YT",Hydrodata$sites)] <- "YT"
Hydrodata$sites[grepl("MOT",Hydrodata$sites)] <- "MOT"
Hydrodata$sites[grepl("NCT",Hydrodata$sites)] <- "NCT"
Hydrodata$sites[grepl("NN",Hydrodata$sites)] <- "NN"
Hydrodata$HydroYear=as.factor(Hydrodata$HydroYear)


mine<-aggregate(Hydrodata,by = list(Hydrodata$sites,Hydrodata$HydroYear,FUN=sum)

Outputv2=Output[which(Hydrodata$wetyear=="FALSE"),]