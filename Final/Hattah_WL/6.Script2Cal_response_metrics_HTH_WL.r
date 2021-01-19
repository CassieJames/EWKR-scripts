##################################################################################################################
#### Script to generate response metrics - ammended 18/11/2019 to include Tda exoitc and natives as a separate metric and remove Tda from amphibious species
library(vegan)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
mydata=read.csv("Spp_site_year_transect matrix HTH_WL_May 2019.csv",row.names=1) # save data out


date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
env.data<-read.csv("Hattah Lakes wetlands transect based env data.csv") # save data out to chec

env.data=env.data[,c("Unique_site_year", "Site.ID.x","Unique_site_year_season","Date.of.collection.x","Elevation", "Easting", "Northing", "Season","Wetland","d90", "d365","TSLW","TSLW_LT","d3Mon_wet",
"d3Mon_meandepth","d3Mon_drylength", "d1yrs_wet", "d1yrs_meandepth","d1yrs_drylength", "d3yrs_wet", "d3yrs_meandepth","d3yrs_drylength","Freq_d1", "Freq_d3", "Freq_d5", "Freq_d10", "Freq_ALL", "WaterYr", 
"Inundated", "d5yearsFF", "d10yearsFF", "d20yearsFF", "d30yearsFF", "MeanTemp90", "MaxTemp90","MinTemp90", "MeanTemp365", "MaxTemp365", "MinTemp365",
"VEG_CLASS", "Max.CTF.30yrs", "P.CTF.30yrs","Max.CTF.20yrs", "P.CTF.20yrs","Max.CTF.10yrs", "P.CTF.10yrs","Max.CTF.5yrs", "P.CTF.5yrs")]


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

Amph=c("ATw","ATl","ARp","ARf","ATe","S","Se")
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
TdagroupNats=Tdagroup[which(Tdagroup$Weed.status == "FALSE"),]
TdagroupExotic=Tdagroup[which(Tdagroup$Weed.status == "TRUE"),]

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
TdagroupNats=TdagroupNats[,-which(names(TdagroupNats) %in% c("Species","Final_PFG_allocation", "Weed.status"))]
TdagroupExotic=TdagroupExotic[,-which(names(TdagroupExotic) %in% c("Species","Final_PFG_allocation", "Weed.status"))]


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
TdagroupNats=as.data.frame(t(TdagroupNats))
TdagroupExotic=as.data.frame(t(TdagroupExotic))

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
Tda_diversity=diversity(TdagroupNats, index = "shannon", MARGIN = 1, base = exp(1))

WetgroupNatsMet=as.data.frame(rowSums(WetgroupNats))
WetgroupExoticMet=as.data.frame(rowSums(WetgroupExotic))
DrygroupNatsMet=as.data.frame(rowSums(DrygroupNats))
DrygroupExoticMet=as.data.frame(rowSums(DrygroupExotic))
TdagroupNatsMet=as.data.frame(rowSums(TdagroupNats))
TdagroupExoticMet=as.data.frame(rowSums(TdagroupExotic))

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
TdagroupNats[TdagroupNats>0] <-1
TdagroupExotic[TdagroupExotic>0] <-1

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
TdagroupNatsRich=as.data.frame(rowSums(TdagroupNats))
TdagroupExoticRich=as.data.frame(rowSums(TdagroupExotic))

ATlgroupRich=as.data.frame(rowSums(ATlgroup))
ATegroupRich=as.data.frame(rowSums(ATegroup))
ARfgroupRich=as.data.frame(rowSums(ARfgroup))
ARpgroupRich=as.data.frame(rowSums(ARpgroup))
ATwgroupRich=as.data.frame(rowSums(ATwgroup))
TdagroupRich=as.data.frame(rowSums(Tdagroup))
TdrgroupRich=as.data.frame(rowSums(Tdrgroup))

HT_WL_metrics <-cbind(WetgroupNatsMet,WetgroupExoticMet,DrygroupNatsMet,DrygroupExoticMet,WetgroupNatsRich,WetgroupExoticRich,DrygroupNatsRich,DrygroupExoticRich,TdagroupExoticRich,TdagroupNatsRich, wet_diversity, dry_diversity,Tda_diversity,
ATlgroupMet,ATegroupMet,ARfgroupMet,ARpgroupMet,ATwgroupMet,TdagroupMet,TdrgroupMet,ATlgroupRich,ATegroupRich,ARfgroupRich,ARpgroupRich,ATwgroupRich,TdagroupRich,TdrgroupRich)

rownames(HT_WL_metrics)=gsub("LHT","LHAT",rownames(HT_WL_metrics))
rownames(HT_WL_metrics)=gsub("CHT","CCS",rownames(HT_WL_metrics))
rownames(HT_WL_metrics)=gsub("CCNT","NCT",rownames(HT_WL_metrics))
rownames(HT_WL_metrics)=gsub("KRT","KT",rownames(HT_WL_metrics))

colnames(HT_WL_metrics)<-c("Wet_Natives", "Wet_Exotics", "Terr_Natives", "Terr_Exotics","WetNatRich", "WetExoticRich","TerrNatRich", "TerrExoticRich" ,"TdaNatRich", "TdaExoticRich", "Wet_diversity", "Dry_diversity","Tda_diversity",
"ATl", "ATe", "ARf", "ARp", "ATw", "Tda", "Tdr","ATlRich", "ATeRich", "ARfRich", "ARpRich", "ATwRich", "TdaRich", "TdrRich" )

rownames(HT_WL_metrics)=gsub("CHT","CCS",rownames(HT_WL_metrics))
rownames(HT_WL_metrics)=gsub("LHT","LHAT",rownames(HT_WL_metrics))

#########################################################################################################

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)




mydata_metrics <- merge(env.data,HT_WL_metrics, by.x="Unique_site_year_season", by.y="row.names")

write.csv(mydata_metrics, file = "Hattah wetlands response by metrics.csv") # save data out to chec

