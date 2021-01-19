##################################################################################################################
#### Script to generate response metrics
library(vegan)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
mydata=read.csv("Spp_site_year_transect matrix HTH_WL_May 2019.csv",row.names=1) # save data out


date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
env.data<-read.csv("Hattah Lakes wetlands transect based env data.csv") # save data out to chec

#########################################################################################################
#### Separate data into broad functional groups and exotic status

species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended
mycodes=mycodes[!duplicated(mycodes[,c('sp_code_simple')]),] # remove duplicates

myspecies=as.data.frame(colnames(mydata))
colnames(myspecies)="Species"
fgrps=merge(myspecies,mycodes,by.x="Species",by.y="sp_code_simple",all.x=TRUE)

fgrps=fgrps[,c("Species","Final_PFG_allocation","Weed.status")]
fgrps$Species=as.character(fgrps$Species)

twetlandALL=as.data.frame(t(mydata))

fgrpswetlands=merge(fgrps,twetlandALL,by.x="Species",by.y="row.names",all.y=TRUE)

Amph=c("ARp","ARf","S","Se","ATl","ATe","ATw")
ATl=c("ATl")
ATe=c("ATe")
ARf=c("ARf")
ARp=c("ARp")
ATw=c("ATw")
Tda=c("Tda")
Tdr=c("Tdr")
Terre=c("Tdr")

ATlgroup=fgrpswetlands[which(fgrpswetlands$Final_PFG_allocation %in% ATl),]
ATegroup=fgrpswetlands[which(fgrpswetlands$Final_PFG_allocation %in% ATe),]
ARfgroup=fgrpswetlands[which(fgrpswetlands$Final_PFG_allocation %in% ARf),]
ARpgroup=fgrpswetlands[which(fgrpswetlands$Final_PFG_allocation %in% ARp),]
ATwgroup=fgrpswetlands[which(fgrpswetlands$Final_PFG_allocation %in% ATw),]
Tdagroup=fgrpswetlands[which(fgrpswetlands$Final_PFG_allocation %in% Tda),]
Tdrgroup=fgrpswetlands[which(fgrpswetlands$Final_PFG_allocation %in% Tdr),]

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
rownames(ATlgroup)=ATlgroup$Species
rownames(ATegroup)=ATegroup$Species
rownames(ARfgroup)=ARfgroup$Species
rownames(ARpgroup)=ARpgroup$Species
rownames(ATwgroup)=ATwgroup$Species
rownames(Tdagroup)=Tdagroup$Species
rownames(Tdrgroup)=Tdrgroup$Species


WetgroupNats=WetgroupNats[,-which(names(WetgroupNats) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
WetgroupExotic=WetgroupExotic[,-which(names(WetgroupExotic) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
DrygroupNats=DrygroupNats[,-which(names(DrygroupNats) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
DrygroupExotic=DrygroupExotic[,-which(names(DrygroupExotic) %in% c("Species","Final_PFG_allocation", "Weed.status"))]

ATlgroup=ATlgroup[,-which(names(ATlgroup) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
ATegroup=ATegroup[,-which(names(ATegroup) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
ARfgroup=ARfgroup[,-which(names(ARfgroup) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
ARpgroup=ARpgroup[,-which(names(ARpgroup) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
ATwgroup=ATwgroup[,-which(names(ATwgroup) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
Tdagroup=Tdagroup[,-which(names(Tdagroup) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
Tdrgroup=Tdrgroup[,-which(names(Tdrgroup) %in% c("Species","Final_PFG_allocation", "Weed.status"))]

WetgroupNats=as.data.frame(t(WetgroupNats))
WetgroupExotic=as.data.frame(t(WetgroupExotic))
DrygroupNats=as.data.frame(t(DrygroupNats))
DrygroupExotic=as.data.frame(t(DrygroupExotic))
ATlgroup=as.data.frame(t(ATlgroup))
ATegroup=as.data.frame(t(ATegroup))
ARfgroup=as.data.frame(t(ARfgroup))
ARpgroup=as.data.frame(t(ARpgroup))
ATwgroup=as.data.frame(t(ATwgroup))
Tdagroup=as.data.frame(t(Tdagroup))
Tdrgroup=as.data.frame(t(Tdrgroup))


TC=rowSums(WetgroupNats)+rowSums(WetgroupExotic)+rowSums(DrygroupNats)+rowSums(DrygroupExotic)


wet_diversity=diversity(WetgroupNats, index = "shannon", MARGIN = 1, base = exp(1))
dry_diversity=diversity(DrygroupNats, index = "shannon", MARGIN = 1, base = exp(1))

WetgroupNatsMet=as.data.frame(rowSums(WetgroupNats))
WetgroupExoticMet=as.data.frame(rowSums(WetgroupExotic))
DrygroupNatsMet=as.data.frame(rowSums(DrygroupNats))
DrygroupExoticMet=as.data.frame(rowSums(DrygroupExotic))
ATlgroupMet=as.data.frame(rowSums(ATlgroup))
ATegroupMet=as.data.frame(rowSums(ATegroup))
ARfgroupMet=as.data.frame(rowSums(ARfgroup))
ARpgroupMet=as.data.frame(rowSums(ARpgroup))
ATwgroupMet=as.data.frame(rowSums(ATwgroup))
TdagroupMet=as.data.frame(rowSums(Tdagroup))
TdrgroupMet=as.data.frame(rowSums(Tdrgroup))


WetgroupNats[WetgroupNats>0] <-1
WetgroupExotic[WetgroupExotic>0] <-1
DrygroupNats[DrygroupNats>0] <-1
DrygroupExotic[DrygroupExotic>0] <-1

ATlgroup[ATlgroup>0] <-1
ATegroup[ATegroup>0] <-1
ARfgroup[ARfgroup>0] <-1
ARpgroup[ARpgroup>0] <-1
ATwgroup[ATwgroup>0] <-1
Tdagroup[Tdagroup>0] <-1
Tdrgroup[Tdrgroup>0] <-1



WetgroupNatsRich=as.data.frame(rowSums(WetgroupNats))
WetgroupExoticRich=as.data.frame(rowSums(WetgroupExotic))
DrygroupNatsRich=as.data.frame(rowSums(DrygroupNats))
DrygroupExoticRich=as.data.frame(rowSums(DrygroupExotic))

ATlgroupRich=as.data.frame(rowSums(ATlgroup))
ATegroupRich=as.data.frame(rowSums(ATegroup))
ARfgroupRich=as.data.frame(rowSums(ARfgroup))
ARpgroupRich=as.data.frame(rowSums(ARpgroup))
ATwgroupRich=as.data.frame(rowSums(ATwgroup))
TdagroupRich=as.data.frame(rowSums(Tdagroup))
TdrgroupRich=as.data.frame(rowSums(Tdrgroup))



HT_WL_metrics <-cbind(WetgroupNatsMet,WetgroupExoticMet,DrygroupNatsMet,DrygroupExoticMet,WetgroupNatsRich,WetgroupExoticRich,DrygroupNatsRich,DrygroupExoticRich, wet_diversity, dry_diversity,
ATlgroupMet,ATegroupMet,ARfgroupMet,ARpgroupMet,ATwgroupMet,TdagroupMet,TdrgroupMet,ATlgroupRich,ATegroupRich,ARfgroupRich,ARpgroupRich,ATwgroupRich,TdagroupRich,TdrgroupRich)


rownames(HT_WL_metrics)=gsub("LHT","LHAT",rownames(HT_WL_metrics))
rownames(HT_WL_metrics)=gsub("CHT","CCS",rownames(HT_WL_metrics))
rownames(HT_WL_metrics)=gsub("CCNT","NCT",rownames(HT_WL_metrics))
rownames(HT_WL_metrics)=gsub("KRT","KT",rownames(HT_WL_metrics))


colnames(HT_WL_metrics)<-c("Aqua_Natives", "Wet_Exotics", "Terr_Natives", "Terr_Exotics","WetNatRich", "WetExoticRich","TerrNatRich", "TerrExoticRich" , "Wet_diversity", "Dry_diversity",
"ATl", "ATe", "ARf", "ARp", "ATw", "Tda", "Tdr","ATlRich", "ATeRich", "ARfRich", "ARpRich", "ATwRich", "TdaRich", "TdrRich" )
#########################################################################################################

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
mydata_metrics <- merge(env.data,HT_WL_metrics, by.x="Unique_site_year_season", by.y="row.names")


write.csv(mydata_metrics, file = "Hattah wetlands response by metrics_aquatics.csv") # save data out to chec

#########################################################################################################
# Method to work out change in vegetation between years - would need to eliminate first year of survey from analysis

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
output=read.csv("Hattah wetlands response by metrics.csv") # load data
output$Wet_Natives_diff <-0

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")

for(w in wetlist) { # 

tdata=output[grep(w,(output$Site.ID.x)),]
sites=unique(tdata$Site.ID.x)

for(s in sites) { # 
s=as.character(s)
ttdata=tdata[grep(s,(tdata$Site.ID.x), fixed=TRUE),]
years=unique(ttdata$WaterYr)

for (i in 1:length(years)){

tttdata=ttdata$Wet_Natives[i+1]-ttdata$Wet_Natives[i]

output[grep(paste(s,"_",years[i]+1,sep=""),output$Unique_site_year_season, fixed=TRUE),ncol(output)] <- tttdata
output=subset(output, output$WaterYr) # just model positive counts

}}}


