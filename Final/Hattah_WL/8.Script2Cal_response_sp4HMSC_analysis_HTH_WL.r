##################################################################################################################
#### Script to generate response metrics - ammended 18/11/2019 to include Tda exoitc and natives as a separate metric and remove Tda from amphibious species
library(vegan)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
mydata=read.csv("Spp_site_year_transect matrix HTH_WL_May 2019.csv",row.names=1) # save data out


date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
env.data<-read.csv("Hattah Lakes wetlands transect based env data January 2021.csv") # save data out to chec


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



#########################################################################################################

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)


mydata_metrics <- merge(env.data,WetgroupNats, by.x="Unique_site_year_season", by.y="row.names")
mydata_metrics <- merge(env.data,DrygroupNats, by.x="Unique_site_year_season", by.y="row.names")

write.csv(mydata_metrics, file = "Hattah wetlands WetgroupNats sp vs enviro.csv") # save data out to chec

