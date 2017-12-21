# Script to turn data into a site by species matrix, determine species metrics and merge with environmental attributes
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 15th December 2017
###########################################################################################
# This section creates a matrix with the environmental attributes tagged onto the end
# load environmental data from various sources - this is messy because the predictor variables were/are in various locations and states

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL_FGcorrections.csv")) # functional group allocations corrected in veg database (some species had two different FG assignments)
mycodes=data.frame(read.csv("HL_WL_sp_codes.csv")) # creation of unique species codes
mydata=merge(mydata,mycodes,by="Scientific.name")
mydata <- within(mydata, Date.of.collection <- as.Date(as.character(Date.of.collection), format = "%d/%m/%Y")) # ensure dates are recognised as dates!

# correct site codes so they match environmental data :(
mydata$Site.ID=gsub("CCNT","NCT",mydata$Site.ID) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie - my attempt to make everything consistent
mydata$Site.ID=gsub("KRT","KT",mydata$Site.ID)
mydata$Unique_site_year=gsub("CCNT","NCT",mydata$Unique_site_year) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie
mydata$Unique_site_year=gsub("KRT","KT",mydata$Unique_site_year)

# Import corrected dates (I corrected the original database for the FP but not the WL dataset - not sure why)
date.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv")) # note that in the original vegetation database the dates of collection were incorrect for a number of years
newdates$Date.of.collection <- as.Date(newdates$Date.of.collection, format="%d/%m/%Y") 
newdates$Unique_site_year=gsub("KRT","KT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CCNT","NCT",newdates$Unique_site_year)
newdates=newdates[,c(1,3,4,5,6)]
mydata2=merge(mydata, newdates, by="Unique_site_year" ) # 

#Hydrology
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd(data.dir)
Hydrodata=data.frame(read.csv("Flood_HTH_WL_pumps_corrected_60 day interval using RIMFIM CTFV2.csv")) # load corrected hydro data with duplicates and date errors removed (V2 is version where contradictions in 'inundation' have been corrected)
Hydrodata$Unique_site_year=gsub("KRT","KT",Hydrodata$Unique_site_year)
Hydrodata$Unique_site_year=gsub("CCNT","NCT",Hydrodata$Unique_site_year)
Hydrodata=Hydrodata[!duplicated(Hydrodata), ] # remove duplicates
Hydrodata=Hydrodata[,c(3,6,8,19,20)] # tidy up and remove extra columns not needed

#Rainfall
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes rainfall/"; setwd(data.dir)
Rainfall.dat=data.frame(read.csv("Rainfall_HTH_WL.csv")) # load rainfall data which is already sorted into dates x rainfall metrics
Rainfall.dat <- within(Rainfall.dat, Date <- as.Date(as.character(Date), format = "%m/%d/%Y")) # ensure dates are recognised in rainfall dat
Rainfall.dat=Rainfall.dat[!duplicated(Rainfall.dat), ] # remove duplicates
Rainfall.dat=Rainfall.dat[,c(2:6)]

#Temperature
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes temperature/"; setwd(data.dir)
Temperature.dat=data.frame(read.csv("Temperature_HTH_WL.csv")) # load temperature data which is already sorted into dates x rainfall metrics
Temperature.dat <- within(Temperature.dat, Date <- as.Date(as.character(Date), format = "%Y-%m-%d")) # ensure dates are recognised in rainfall dat
Temperature.dat=Temperature.dat[,c(2:14)]
Temperature.dat=Temperature.dat[!duplicated(Temperature.dat), ] # remove duplicates

#Veg overstorey structure and location
data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes veg structure/"; setwd(data.dir)
VegStructure.dat=data.frame(read.csv("Hattah WL Veg Structure.csv")) # load veg structure data
VegStructure.dat<-VegStructure.dat[,c(2,14)] 

#Merge various datasets
mydata_env=merge(mydata2, Rainfall.dat, by.x="Date.of.collection.y", by.y="Date", all.x=TRUE) # merge species data with rainfall data by date to create a unique site-sample date by rainfall metrics table
mydata_env=merge(mydata_env, Hydrodata, by="Unique_site_year", all.x=TRUE) # I am not sure why but when I merge its duplicating selected lines of the database....
mydata_env=merge(mydata_env, Temperature.dat, by.x="Date.of.collection.y", by.y="Date",all.x=TRUE) # merge 
mydata_env=merge(mydata_env, VegStructure.dat, by="Site.ID", all.x=TRUE) # merge  

# Had to change some of the site codes as the Hattah lakes site codes are a subset of some of the Chalka creek and Little Hattah site codes so the matrix was filled in correctly
mydata_env$Unique_site_year=gsub("CHT","CCS",mydata_env$Unique_site_year)
mydata_env$Site.ID=gsub("CHT","CCS",mydata_env$Site.ID)
mydata_env$Unique_site_year=gsub("LHT","LHAT",mydata_env$Unique_site_year)
mydata_env$Site.ID=gsub("LHT","LHAT",mydata_env$Site.ID)
mydata_env_OI=mydata_env[,c("Unique_site_year", "Inundated","d30", "d365", "Flood_frequency", "MaxTemp365", "MinTemp365","VEG_CLASS")] # subset to useful env data
mydata_env_OI=mydata_env_OI[!duplicated(mydata_env_OI), ] # remove duplicates (I think I solved this so probably are not any more duplicates)

# Create species x site matrix and tag on environmental data
specieslist=unique(mydata_env$sp_code)
sitelist=unique(mydata_env$Unique_site_year)
takeout=c("Inundate", "Lea.litt", "Bar.grou") # remove non species
specieslist=specieslist[!specieslist %in% takeout]
Output= matrix(NA,nrow=length(unique(mydata2$Unique_site_year)), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=sort(specieslist)
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata_env[which(mydata_env$Unique_site_year==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]

if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$sp_code==spp)])
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix


Spp_Env_HTH_WL=merge(Output,mydata_env_OI ,by.x="row.names", by.y="Unique_site_year") # merge 

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
write.csv(Spp_Env_HTH_WL , file = "Spp_Env_HTH_WL_Dec 2017.csv") # save data out

###########################################################################################
## Script to aggregate data to wetland and year

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")
yoi=c("_08", "_09","_10","_11", "_12", "_13", "_14", "_16")
tt = expand.grid(wetlist,yoi); tt = paste(tt[,1],tt[,2],sep='_')

tada = matrix(NA,nrow=length(tt),ncol=ncol(Output))#define the output matrix
rownames(tada)=tt
colnames(tada)=colnames(Output)

for(w in wetlist) { # 
tdata=Output[grep(w,rownames(Output)),]

for (yy in yoi) {
ttdata=tdata[grep(yy,rownames(tdata)),]
out=colSums(ttdata)
tada[grep(paste(w,yy,sep="_"),rownames(tada)),] = out
}

}

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_matrix summarised to wetland_HTH_WL_Dec 2017.csv") # save data out















###########################################################################################
#Scripts to calculate various vegetation metrics

######## Exotic metrics ########
# Determine proportion of occurrences that are exotic species

Outputt=as.data.frame(Output)
Outputt=t(Output) # transpose matrix as my brain works better in this direction!
Outputt=as.data.frame(Outputt)
Outputt$Exotic <- NA # create an empty column for the exotic status to fill

for(spp in specieslist) { # loops through species list
exotic.status=unique(mydata[mydata$sp_code==spp,c("Weed.status")]) # extract weed status from original mydata frame
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

# rem the above code has also summed the functional group designation columns!

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


####### H index and species richness ########

library(vegan)

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
Spp_Env_matrix_HTH_WL=read.csv("Spp_Env_HTH_WL_Dec 2017.csv") # read data
spp.matrix=Spp_Env_matrix_HTH_WL[,c(2:239)] # this is lazy coding so need to be careful it picks up all the species
rownames(spp.matrix)=spp.matrix$Row.names
spp.matrix=spp.matrix[,-c(1)]
spp.matrix.pa=spp.matrix
spp.matrix.pa[ spp.matrix.pa> 0]<- 1 # make a copy of the spp matrix and turn into p/a data
Richness=as.data.frame(rowSums(spp.matrix.pa))
colnames(Richness)=c("Richness")
H.index=(as.data.frame(diversity(spp.matrix, index="shannon")))
colnames(H.index)=c("H_index")

###### Merge metrics with environmental data and save out for analysis

Final_Hattah_WL=merge(tdata, ttdata, by.x="row.names", by.y="row.names", all.x=TRUE) #
Final_Hattah_WL=merge(Final_Hattah_WL, H.index, by.x="Row.names", by.y="row.names", all.x=TRUE) #
Final_Hattah_WL=merge(Final_Hattah_WL, Richness, by.x="Row.names", by.y="row.names", all.x=TRUE) #
Final_Hattah_WL=merge(Final_Hattah_WL, mydata_env_OI, by.x="Row.names", by.y="Unique_site_year", all.x=TRUE) # merge with restricted enviro data with duplicates removed
write.csv(Final_Hattah_WL , file = "Final_Metrics_HTH_WL.csv") # save data out

###### Some basic plots of the responses

envdata=read.csv("Final_Metrics_HTH_WL.csv") 
envdata$Terrestrial=envdata$T+envdata$Tda+envdata$Tdr # Here I have just grouped the functional groups into coarser categories
envdata$Aquatic=envdata$A+envdata$Atw+envdata$Atl+envdata$Ate+envdata$Arp+envdata$F+envdata$S
envdata$Aquatic_prop =(envdata$Aquatic/(envdata$Terrestrial+envdata$Aquatic))*100
envdata$Terrestrial_prop =(envdata$Terrestrial/(envdata$Terrestrial+envdata$Aquatic))*100
envdata$Total_abund=envdata$Terrestrial+envdata$Aquatic

library(ggplot2)
library(gridExtra)
image.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Plots/"

png(paste(image.dir,'Response Metrics hists WL.png',sep=''), width=1500, height=1500, units="px", res=200)

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
