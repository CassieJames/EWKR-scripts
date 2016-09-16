# Script to turn data into a site by species matrix, determine species metrics and merge with environmental attributes
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016

# Import data
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))
mycodes=data.frame(read.csv("HL_FP_sp_codes.csv"))
mydata=merge(mydata,mycodes,by.x="Scientific_name", by.y="Scientific.name")


###########################################################################################
# This section creates a matrix with the environmental attributes tagged onto the end
# load environmental data from various sources - this is messy because the predictor variables were/are in various locations and states

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))
mycodes=data.frame(read.csv("HL_FP_sp_codes.csv"))
mydata=merge(mydata,mycodes,by.x="Scientific_name", by.y="Scientific.name")

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/"; setwd(data.dir)
Hydrodata=data.frame(read.csv("Flood_HTH_FP_pumps_corrected.csv")) # load corrected hydro data with duplicates and date errors removed
Rainfall.dat=data.frame(read.csv("Rainfall_HTH_FP.csv")) # load rainfall data which is already sorted into dates x rainfall metrics
Rainfall.dat <- within(Rainfall.dat, Date <- as.Date(as.character(Date), format = "%d/%m/%Y")) # ensure dates are recognised in rainfall dat
GIS.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/GIS stuff/"; setwd(GIS.dir)
HTH_FP_locs=data.frame(read.csv("HTH_FP_locations.csv")) # load site longs and lats
Hydrodata=Hydrodata[!duplicated(Hydrodata), ] # remove duplicates
mydata <- within(mydata, Date.of.collection <- as.Date(as.character(Date.of.collection), format = "%d/%m/%Y")) # ensure dates are still recognised!
Rainfall.dat=Rainfall.dat[!duplicated(Rainfall.dat), ] # remove duplicates
mydata_env=merge(mydata, Rainfall.dat, by.x="Date.of.collection", by.y="Date") # merge species data with rainfall data by date to create a unique site-sample date by rainfall metrics table
mydata_rainfall=mydata_env[,c("Unique_site_year", "d30", "d90", "d180", "d365")]
mydata_rainfall=unique(mydata_rainfall)
mydata_env=merge(mydata_rainfall, Hydrodata, by="Unique_site_year") # merge with hydrodata - I have used merge as this is safer in case the rows are not in the same order
mydata_env=merge(mydata_env, HTH_FP_locs, by="Site.ID") # merge with location info
wrc.data=mydata[,c("Site.ID", "WRC")]# extracts WRC and attributes to site
wrc.data=wrc.data[!duplicated(wrc.data), ]
mydata_env=merge(mydata_env, wrc.data, by="Site.ID") # merge with location info

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes temperature/"; setwd(data.dir)
Temperature.dat=data.frame(read.csv("Temperature_HTH_FP.csv")) # load temperature data which is already sorted into dates x rainfall metrics
Temperature.dat <- within(Temperature.dat, Date <- as.Date(as.character(Date), format = "%Y-%m-%d")) # ensure dates are recognised in rainfall dat
Temperature.dat=Temperature.dat[,c(2:14)]
Temperature.dat=Temperature.dat[!duplicated(Temperature.dat), ] # remove duplicates
mydata_env <- within(mydata_env, Date.of.collection <- as.Date(as.character(Date.of.collection), format = "%d/%m/%Y")) # ensure dates are still recognised!
mydata_env=merge(mydata_env, Temperature.dat, by.x="Date.of.collection", by.y="Date") # merge 
 
# Create species x site matrix and tag on environmental data
specieslist=unique(mydata$sp_code)
sitelist=unique(mydata$Unique_site_year)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
specieslist=specieslist[!specieslist %in% takeout]
Output= matrix(NA,nrow=length(unique(mydata$Unique_site_year)), ncol=length(specieslist))
rownames(Output)=unique(mydata$Unique_site_year)
colnames(Output)=sort(specieslist)
Output=as.data.frame(Output)
Output=merge(Output, mydata_env, by.x="row.names", by.y="Unique_site_year", all.x=TRUE) # add rainfall metrics to empty site x species matrix

for(s in sitelist) { # Fill matrix
tdata=mydata[which(mydata$Unique_site_year==s),]
sitesp=unique(tdata$sp_code)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$sp_code==spp)])
Output[grep(s, Output$Row.names),grep(spp,colnames(Output))] <- abund
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix


Spp_Env_matrix_HTH_FP=Output # take a copy
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
write.csv(Spp_Env_matrix_HTH_FP , file = "Spp_Env_matrix_HTH_FP.csv") # save data out

###########################################################################################
#Scripts to calculate various vegetation metrics

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))
mycodes=data.frame(read.csv("HL_FP_sp_codes.csv"))
mydata=merge(mydata,mycodes,by.x="Scientific_name", by.y="Scientific.name")

#create matrix for results (species x sites)
specieslist=unique(mydata$sp_code)
sitelist=unique(mydata$Unique_site_year)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
specieslist=specieslist[!specieslist %in% takeout]
sitelist=unique(mydata$Unique_site_year)
Output= matrix(NA,nrow=length(unique(mydata$Unique_site_year)), ncol=length(specieslist))
rownames(Output)=unique(mydata$Unique_site_year)
colnames(Output)=sort(specieslist)

for(s in sitelist) {
tdata=mydata[which(mydata$Unique_site_year==s),]
sitesp=unique(tdata$sp_code)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$sp_code==spp)])
Output[grep(s, rownames(Output)),grep(spp,colnames(Output))] <- abund
}
}
Output[is.na(Output)] <- 0

######## Exotic metrics ########
# Determine proportion of occurrences that are exotic species

Outputt=as.data.frame(Output)
Outputt=t(Output) # transpose matrix as my brain works better in this direction!
Outputt=as.data.frame(Outputt)
Outputt$Exotic <- NA # create an empty column for the exotic status to fill

for(spp in specieslist) { # loops through species list
exotic.status=unique(mydata[mydata$sp_code==spp,c("Weed")]) # extract weed status from original mydata frame
Outputt[agrep(spp, rownames(Outputt)),"Exotic"] <- exotic.status} # agrep gives a lazy match - it seems to work where the exact match was falling over for a couple of species - although why I don't know :(

exotic=colSums(Outputt[which(Outputt$Exotic==2),]) # sum all rows for each column for which the species is exotic
native=colSums(Outputt[which(Outputt$Exotic==1),]) # sum all rows for each column for which the species is native
total_abund=colSums(Outputt)
exotic_prop=(exotic/total_abund)*100
tdata=rbind(exotic,native)
tdata=rbind(tdata,total_abund)
tdata=rbind(tdata,exotic_prop)

######## Functional Group ########
# Determine abundances of different functional groups

Outputt=as.data.frame(Output)
Outputt=t(Output) # transpose matrix as my brain works better in this direction!
Outputt=as.data.frame(Outputt)
Outputt$Tdr <- NA # create empty columns
Outputt$Tda <- NA
Outputt$T <- NA
Outputt$A <- NA
Outputt$Atw <- NA
Outputt$Atl <- NA
Outputt$Ate <- NA
Outputt$Arp <- NA
Outputt$F <- NA
Outputt$S <- NA
specieslist=colnames(Output)

FGs=c("T", "Tda", "Tdr", "A", "Arp","Atw", "Atl", "Ate", "Arp", "F", "S") # identify functional groups of interest - remove leaf litter and bare ground
for(spp in specieslist) { # loops through species list
Fung.group=unique(mydata[mydata$sp_code==spp,c("Functional.group")]) # for each species identify its functional group
if(Fung.group %in% FGs){Outputt[agrep(spp, rownames(Outputt)),match(Fung.group, colnames(Outputt))] <- 1} # is the functional group of the species is in FGs list then identify the species and match the correct functional group
} # agrep gives a lazy match - it seems to work where the exact match was falling over for a couple of species

T=colSums(Outputt[which(Outputt$T==1),]) # Sum all rows for which the functional group equals T and so on...
Tdr=colSums(Outputt[which(Outputt$Tdr==1),])
Tda=colSums(Outputt[which(Outputt$Tda==1),])
A=colSums(Outputt[which(Outputt$A==1),])
Atw=colSums(Outputt[which(Outputt$Atw==1),])
Atl=colSums(Outputt[which(Outputt$Atl==1),])
Ate=colSums(Outputt[which(Outputt$Ate==1),])
Arp=colSums(Outputt[which(Outputt$Arp==1),])
F=colSums(Outputt[which(Outputt$F==1),])
S=colSums(Outputt[which(Outputt$S==1),])


ttdata=rbind(T,Tdr) # lazy binding!
ttdata=rbind(ttdata,Tda)
ttdata=rbind(ttdata,A)
ttdata=rbind(ttdata,Atw)
ttdata=rbind(ttdata,Atl)
ttdata=rbind(ttdata,Ate)
ttdata=rbind(ttdata,Arp)
ttdata=rbind(ttdata,F)
ttdata=rbind(ttdata,S)

tdata=t(tdata) # transpose matrix
REMOVE=c("Exotic")
tdata=tdata[!rownames(tdata) %in% REMOVE, ] # remove exotic row from bottom of data frame
ttdata=t(ttdata)
ttdata=ttdata[!rownames(ttdata) %in% FGs, ] # remove functional rows from bottom of data frame
Final_metrics=merge(tdata,ttdata,by='row.names')

######## H index and species richness ########
# Diversity metric calculated to do trial analysis as it has quite a simple distribution (gaussian)
library(vegan)

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
Spp_Env_matrix_HTH_FP=read.csv("Spp_Env_matrix_HTH_FP.csv") # read data
spp.matrix=Spp_Env_matrix_HTH_FP[,c(2:267)]
rownames(spp.matrix)=spp.matrix$Row.names
spp.matrix=spp.matrix[,-c(1)]
spp.matrix.pa=spp.matrix
spp.matrix.pa[ spp.matrix.pa> 0]<- 1 # make a copy of the spp matrix and turn into p/a data
Richness=as.data.frame(rowSums(spp.matrix.pa))
colnames(Richness)=c("Richness")
H.index=(as.data.frame(diversity(spp.matrix, index="shannon")))
colnames(H.index)=c("H_index")

###### Merge metrics with environmental data and save out for analysis

Final_Hattah_FP=merge(Final_metrics, H.index, by.x="Row.names", by.y="row.names", all.x=TRUE) #
Final_Hattah_FP=merge(Final_Hattah_FP, Richness, by.x="Row.names", by.y="row.names", all.x=TRUE) #
Final_Hattah_FP=merge(Final_Hattah_FP, mydata_env, by.x="Row.names", by.y="Unique_site_year", all.x=TRUE) #
write.csv(Final_Hattah_FP , file = "Final_Metrics_Hattah_FP_with rich.csv") # save data out

###### Some basic plots of the responses

envdata=read.csv("Final_Metrics_Hattah_FP_with rich.csv") # load data with corrections to diversity, richness and abundance
envdata$Terrestrial=envdata$T+envdata$Tda+envdata$Tdr
envdata$Aquatic=envdata$A+envdata$Atw+envdata$Atl+envdata$Ate+envdata$Arp+envdata$F+envdata$S
envdata$Aquatic_prop =(envdata$Aquatic/(envdata$Terrestrial+envdata$Aquatic))*100
envdata$Terrestrial_prop =(envdata$Terrestrial/(envdata$Terrestrial+envdata$Aquatic))*100
envdata$Total_abund=envdata$Terrestrial+envdata$Aquatic


png(paste(image.dir,'Response Metrics hists.png',sep=''), width=1500, height=1500, units="px", res=200)

prich <- ggplot(as.data.frame(envdata$Richness), aes(x=envdata$Richness),xlab="Richness")+ geom_histogram(bins=20)
prich <- prich + geom_histogram(bins=20) + labs(x="Richness", y="Count")

pH <- ggplot(as.data.frame(envdata$H_index), aes(x=envdata$H_index))
pH <- pH + geom_histogram(bins=20) + labs(x="H index", y="Count")

pabund <- ggplot(as.data.frame(envdata$total_abund), aes(x=envdata$total_abund))
pabund <- pabund + geom_histogram(bins=20) + labs(x="Total count", y="Count")

pexotic <- ggplot(as.data.frame(envdata$exotic), aes(x=envdata$exotic))
pexotic <- pexotic + geom_histogram(bins=20) + labs(x="Exotic count", y="Count")

pexotic_prop <- ggplot(as.data.frame(envdata$exotic_prop), aes(x=envdata$exotic_prop))
pexotic_prop <- pexotic_prop + geom_histogram(bins=20) + labs(x="Exotic %", y="Count")

pAquatic_prop <- ggplot(as.data.frame(envdata$Aquatic_prop), aes(x=envdata$Aquatic_prop))
pAquatic_prop <- pAquatic_prop + geom_histogram(bins=20) + labs(x="Aquatic %", y="Count")

pTerrestrial_prop <- ggplot(as.data.frame(envdata$Terrestrial_prop), aes(x=envdata$Terrestrial_prop))
pTerrestrial_prop <- pTerrestrial_prop + geom_histogram(bins=20) + labs(x="Terrestrial %", y="Count")

grid.arrange(prich, pH,pabund,pexotic,pAquatic_prop,pTerrestrial_prop, ncol = 2)


dev.off()
