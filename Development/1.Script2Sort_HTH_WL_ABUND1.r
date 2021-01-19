# Script to turn data into a site by species matrix, determine species metrics and merge with environmental attributes
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 15th December 2017
###########################################################################################
# This section creates a matrix with the environmental attributes tagged onto the end
# load environmental data from various sources - this is messy because the predictor variables were/are in various locations and states

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL_FGcorrections.csv")) # functional group allocations corrected in veg database (some species had two different FG assignments)
species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include all sites

mydata=merge(mydata,mycodes,by="Scientific.name",all.x=TRUE)
#write.csv(mydata , file = "Species code check matrix HTH_WL_August 2018.csv") # save data out to check list
mydata <- within(mydata, Date.of.collection <- as.Date(as.character(Date.of.collection), format = "%d/%m/%Y")) # ensure dates are recognised as dates!

# correct site codes so they match environmental data :(
mydata$Site.ID=gsub("CCNT","NCT",mydata$Site.ID) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie - my attempt to make everything consistent
mydata$Site.ID=gsub("KRT","KT",mydata$Site.ID)
mydata$Unique_site_year=gsub("CCNT","NCT",mydata$Unique_site_year) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie
mydata$Unique_site_year=gsub("KRT","KT",mydata$Unique_site_year)

# Import corrected dates (I corrected the original database for the FP but not the WL dataset - not sure why)
date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 

newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv")) # note that in the original vegetation database the dates of collection were incorrect for a number of years
newdates$Date.of.collection <- as.Date(newdates$Date.of.collection, format="%d/%m/%Y") 
newdates$Unique_site_year=gsub("KRT","KT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CCNT","NCT",newdates$Unique_site_year)
newdates=newdates[,c(1,3,4,5,6)]
mydata2=merge(mydata, newdates, by="Unique_site_year" ) # 

#Assign season to dates

library(zoo)
      yq <- as.yearqtr(as.yearmon(mydata2$Date.of.collection.y , "%m/%d/%Y") + 1/12)
      mydata2$Season <- factor(format(yq, "%q"), levels = 1:4, 
      labels = c("SU", "AU", "WI", "SP"))

mydata2$Unique_site_year_season <-paste(mydata2$Unique_site_year,"_",mydata2$Season, sep="")
	  

# Create species x site matrix 
specieslist=unique(mydata2$sp_code_simple) 
specieslist <- specieslist[!is.na(specieslist)]
sitelist=unique(mydata2$Unique_site_year_season)
takeout=c("Inundate", "Lea.litt", "Bar.grou","NA") # remove non species
specieslist=specieslist[!specieslist %in% takeout]
Output= matrix(NA,nrow=length(unique(mydata2$Unique_site_year_season)), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=sort(specieslist)
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata2[which(mydata2$Unique_site_year_season==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code_simple)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]
sitesp <- sitesp[!is.na(sitesp)]

if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$sp_code_simple==spp)])
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 

rownames(Output)=paste("HAT_",rownames(Output),sep="")
rownames(Output)=gsub("_13","_2013",rownames(Output))
rownames(Output)=gsub("_16","_2016",rownames(Output))
rownames(Output)=gsub("_15","_2015",rownames(Output))
rownames(Output)=gsub("_14","_2014",rownames(Output))
rownames(Output)=gsub("_12","_2012",rownames(Output))
rownames(Output)=gsub("_11","_2011",rownames(Output))
rownames(Output)=gsub("_10","_2010",rownames(Output))
rownames(Output)=gsub("_09","_2009",rownames(Output))
rownames(Output)=gsub("_08","_2008",rownames(Output))
rownames(Output)=gsub("_07","_2007",rownames(Output))
rownames(Output)=gsub("_06","_2006",rownames(Output))


write.csv(Output , file = "Spp_site_year_transect matrix HTH_WL_May 2019.csv") # save data out

########################################################################################### 
#### Script to summarise data for each wetland

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")
yoi=c("2000","2001","2002","2003","2004","2005","2006","2007","2008", "2009","2010","2011", "2012", "2013", "2014", "2015","2016")
season=c("SU", "AU","WI","SP")

tt = expand.grid(wetlist); tt = paste("HAT_",tt[,1],sep="")

tada = matrix(NA,nrow=length(tt),ncol=24)#define the output matrix
rownames(tada)=tt
colnames(tada)=c("Totaln","Transectn","Elevn","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","SU","AU","WI","SP")

for(w in wetlist) { # 

tdata=Output[grep(w,rownames(Output)),]
Totaln=nrow(tdata)
UniTs =length(unique(str_extract(rownames(tdata), "\\d")))
Elevs=vapply(strsplit(rownames(tdata),"\\+"), `[`, 2, FUN.VALUE=character(1))
Elevsminus=vapply(strsplit(rownames(tdata),"\\-"), `[`, 2, FUN.VALUE=character(1))
Elevs2=unique(vapply(strsplit(Elevs,"\\_"), `[`, 1, FUN.VALUE=character(1))) # extract number of different elevations for site
Elevs2=length(Elevs2[!is.na(Elevs2)])
Elevsminus2=unique(vapply(strsplit(Elevsminus,"\\_"), `[`, 1, FUN.VALUE=character(1))) # extract number of different elevations for site that are negative
Elevsminus2=length(Elevsminus2[!is.na(Elevsminus2)])
ElevALL=Elevs2+Elevsminus2

tada[grep(w,rownames(tada)),grep("Totaln",colnames(tada))] = Totaln
tada[grep(w,rownames(tada)),grep("Transectn",colnames(tada))] = UniTs
tada[grep(w,rownames(tada)),grep("Elevn",colnames(tada))] = ElevALL

for (yy in yoi) {
yearno=(tdata[grep(yy,rownames(tdata)),])
if (nrow(yearno)>0){
yearno=as.data.frame(yearno)
yearno=yearno[rowSums(yearno!= 0) > 0,]	# remove  rows with no records
}
yearno=nrow(yearno)
tada[grep(w,rownames(tada)),grep(yy,colnames(tada))] = yearno
}

for (seas in season){
seasno=nrow(tdata[grep(seas,rownames(tdata)),])
tada[grep(w,rownames(tada)),grep(seas,colnames(tada))] = seasno
}
}

tada=tada[,colSums(tada!= 0) > 0]	# remove columns with no records

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Data summary table HTH_WL_May 2019.csv") # save data out

########################################################################################### IMPORTANT THIS wont work UNTIL I HAVE SORTED OUT ENVIRONMENTAL DATA ROW HEADINGS
Spp_Env_HTH_WL=merge(Output,mydata_env_OI ,by.x="row.names", by.y="Unique_site_year_season") # merge 

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
write.csv(Spp_Env_HTH_WL , file = "Spp_Env_HTH_WL_May 2019.csv") # save data out

###########################################################################################
## Script to aggregate data to wetland and year

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")
yoi=c("_2008", "_2009","_2010","_2011", "_2012", "_2013", "_2014", "_2016")
tt = expand.grid(wetlist,yoi); tt = paste("HAT_",tt[,1],tt[,2],sep="")

tada = matrix(NA,nrow=length(tt),ncol=ncol(Output))#define the output matrix
rownames(tada)=tt
colnames(tada)=colnames(Output)

for(w in wetlist) { # 
tdata=Output[grep(w,rownames(Output)),]

for (yy in yoi) {
ttdata=tdata[grep(yy,rownames(tdata)),]
maxabund=nrow(ttdata)*15
out=colSums(ttdata)
tada[grep(paste(w,yy,sep=""),rownames(tada)),] = out
}}

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)
tada<-tada[complete.cases(tada), ]

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Spp_site_year matrix HTH_WL_May 2019.csv") # save data out

###########################################################################################
# Summarise Hattah to site

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")

tada2 = matrix(NA,nrow=length(wetlist),ncol=ncol(tada))#define the output matrix
rownames(tada2)=wetlist
colnames(tada2)=colnames(tada)

for(w in wetlist) { # 
tdata=tada [grep(w,rownames(tada)),]
tdata=colSums(tdata)
tada2[grep(w,rownames(tada2)),]<- tdata
}
tada2[is.na(tada2)]<-0

rownames(tada2)=paste("HAT_",rownames(tada2),sep="")

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
write.csv(tada2 , file = "Spp_site matrix HTH_WL_May 2019.csv") # save data out

###########################################################################################
## Script to aggregate data to wetland and year and season

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")
yoi=c("_2008", "_2009","_2010","_2011", "_2012", "_2013", "_2014", "_2016")
season=c("_SU", "_AU","_WI","_SP")

tt = expand.grid(wetlist,yoi,season); tt = paste("HAT_",tt[,1],tt[,2],tt[,3],sep="")

tada = matrix(NA,nrow=length(tt),ncol=ncol(Output))#define the output matrix
rownames(tada)=tt
colnames(tada)=colnames(Output)

for(w in wetlist) { # 
tdata=Output[grep(w,rownames(Output)),]

for (yy in yoi) {
ttdata=tdata[grep(yy,rownames(tdata)),]

for (seas in season){
tttdata=ttdata[grep(seas,rownames(ttdata)),]
maxabund=nrow(tttdata)*15
out=colSums(tttdata)
tada[grep(paste(w,yy,seas,sep=""),rownames(tada)),] = out
}}}

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)
tada<-tada[complete.cases(tada), ]


data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Spp_site_year_season matrix HTH_WL_May 2019.csv") # save data out

###################################################################################################################################
# May 2019 - this creates the dataframe to use for the common species ONLY beta diversity analysis
###################################################################################################################################

Output=ceiling(tada) # first round up as these numbers represent % cover - any detection is a detection
Output[Output > 0] <- 1# Second change to presence/absence because we are just interested in repeated observations over time

HATsites=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")

tada = matrix(NA,nrow=length(HATsites),ncol=ncol(Output))#define the output matrix
rownames(tada)=Hatsites
colnames(tada)=colnames(Output)

for(s in HATsites) { # 
tdata=Output[grep(s,rownames(Output)),]
out=colSums(tdata)
tada[grep(s,rownames(tada)),] <- out
}

rownames(tada)=paste("HAT","_",rownames(tada),sep="")

tada[tada<3] <- 0 # if species recorded less than 5 times turn value to zero

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_matrix Hattah_WL_May 2019 COMMON ONLY.csv") # save data out


################################################################################################################################
# for beta analysis looking at relationships with hydrology
###########################################################################################
## Script to aggregate data to wetland and year and season

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL_FGcorrections.csv")) # functional group allocations corrected in veg database (some species had two different FG assignments)
species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include all sites

mydata=merge(mydata,mycodes,by="Scientific.name",all.x=TRUE)
#write.csv(mydata , file = "Species code check matrix HTH_WL_August 2018.csv") # save data out to check list
mydata <- within(mydata, Date.of.collection <- as.Date(as.character(Date.of.collection), format = "%d/%m/%Y")) # ensure dates are recognised as dates!

# correct site codes so they match environmental data :(
mydata$Site.ID=gsub("CCNT","NCT",mydata$Site.ID) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie - my attempt to make everything consistent
mydata$Site.ID=gsub("KRT","KT",mydata$Site.ID)
mydata$Unique_site_year=gsub("CCNT","NCT",mydata$Unique_site_year) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie
mydata$Unique_site_year=gsub("KRT","KT",mydata$Unique_site_year)

# Import corrected dates (I corrected the original database for the FP but not the WL dataset - not sure why)
date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 

newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv")) # note that in the original vegetation database the dates of collection were incorrect for a number of years
newdates$Date.of.collection <- as.Date(newdates$Date.of.collection, format="%d/%m/%Y") 
newdates$Unique_site_year=gsub("KRT","KT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CCNT","NCT",newdates$Unique_site_year)

newdates=newdates[,c(1,3,4,5,6)]
mydata2=merge(mydata, newdates, by="Unique_site_year" ) # 
mydata2$Unique_site_year=gsub("CHT","CCS",mydata2$Unique_site_year)
mydata2$Unique_site_year=gsub("LHT","LHAT",mydata2$Unique_site_year)

#Assign season to dates

library(zoo)
      yq <- as.yearqtr(as.yearmon(mydata2$Date.of.collection.y , "%m/%d/%Y") + 1/12)
      mydata2$Season <- factor(format(yq, "%q"), levels = 1:4, 
      labels = c("SU", "AU", "WI", "SP"))

mydata2$Unique_site_year_season <-paste(mydata2$Unique_site_year,"_",mydata2$Season, sep="")
	  

# Create species x site matrix 
specieslist=unique(mydata2$sp_code_simple) 
specieslist <- specieslist[!is.na(specieslist)]
sitelist=unique(mydata2$Unique_site_year_season)
takeout=c("Inundate", "Lea.litt", "Bar.grou","NA") # remove non species
specieslist=specieslist[!specieslist %in% takeout]
Output= matrix(NA,nrow=length(unique(mydata2$Unique_site_year_season)), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=sort(specieslist)
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata2[which(mydata2$Unique_site_year_season==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code_simple)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]
sitesp <- sitesp[!is.na(sitesp)]

if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$sp_code_simple==spp)])
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 

rownames(Output)=paste("HAT_",rownames(Output),sep="")
rownames(Output)=gsub("_13","_2013",rownames(Output))
rownames(Output)=gsub("_16","_2016",rownames(Output))
rownames(Output)=gsub("_15","_2015",rownames(Output))
rownames(Output)=gsub("_14","_2014",rownames(Output))
rownames(Output)=gsub("_12","_2012",rownames(Output))
rownames(Output)=gsub("_11","_2011",rownames(Output))
rownames(Output)=gsub("_10","_2010",rownames(Output))
rownames(Output)=gsub("_09","_2009",rownames(Output))
rownames(Output)=gsub("_08","_2008",rownames(Output))
rownames(Output)=gsub("_07","_2007",rownames(Output))
rownames(Output)=gsub("_06","_2006",rownames(Output))



write.csv(Output , file = "Spp_site_year_season matrix HTH_WL_May 2019 for beta analysis.csv") # save data out
