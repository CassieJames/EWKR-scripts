##################################################################################################################
#### Script to generate response metrics

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
mydata=read.csv("Spp_site_year_transect matrix HTH_WL_July 2018.csv",row.names=1) # save data out


date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
env.data<-read.csv("Hattah Lakes wetlands transect based env data.csv") # save data out to chec
env.data=env.data[,c(2,3,6,7,8,9,11,13,14,17,18,19,20,21,22,23,24,25,30,31,38)]

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

TC=rowSums(WetgroupNats)+rowSums(WetgroupExotic)+rowSums(DrygroupNats)+rowSums(DrygroupExotic)

WetgroupNatsMet=as.data.frame(rowSums(WetgroupNats))
WetgroupExoticMet=as.data.frame(rowSums(WetgroupExotic))
DrygroupNatsMet=as.data.frame(rowSums(DrygroupNats))
DrygroupExoticMet=as.data.frame(rowSums(DrygroupExotic))

HT_WL_metrics <-cbind(WetgroupNatsMet,WetgroupExoticMet,DrygroupNatsMet,DrygroupExoticMet)
colnames(HT_WL_metrics)<-c("Wet_Natives", "Wet_Exotics", "Terr_Natives", "Terr_Exotics")
#########################################################################################################

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
mydata_metrics <- merge(env.data,HT_WL_metrics, by.x="Unique_site_year_season", by.y="row.names")
write.csv(mydata_metrics, file = "Hattah wetlands response by metrics.csv") # save data out to chec
