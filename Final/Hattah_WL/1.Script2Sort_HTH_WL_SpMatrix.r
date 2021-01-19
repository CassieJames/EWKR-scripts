###################################################################################################################################
### Script to turn data into a site by species matrix, determine species metrics and merge with environmental attributes
### Written by  C.S.James 
### GNU General Public License .. feel free to use / distribute ... no warranties
### 15th December 2017 updated 27th August 2019
###################################################################################################################################
# Notes: This section creates a matrix with the environmental attributes tagged onto the end
# load environmental data from various sources - this is messy because the predictor variables were/are in various locations and states

# step 1 - load data and correct dates
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL_FGcorrections.csv")) # functional group allocations corrected in veg database (some species had two different FG assignments)
species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include all sites

mydata=merge(mydata,mycodes,by="Scientific.name",all.x=TRUE)
#write.csv(mydata , file = "Species code check matrix HTH_WL_August 2018.csv") # save data out to check list
mydata <- within(mydata, Date.of.collection <- as.Date(as.character(Date.of.collection), format = "%d/%m/%Y")) # ensure dates are recognised as dates!

# correct site codes so they match environmental data :(
mydata$Site.ID=gsub("CCNT","NCT",mydata$Site.ID) 
mydata$Site.ID=gsub("KRT","KT",mydata$Site.ID)
mydata$Site.ID=gsub("CHT","CCS",mydata$Site.ID)
mydata$Site.ID=gsub("LHT","LHAT",mydata$Site.ID)
mydata$Unique_site_year=gsub("CCNT","NCT",mydata$Unique_site_year) 
mydata$Unique_site_year=gsub("KRT","KT",mydata$Unique_site_year)
mydata$Unique_site_year=gsub("CHT","CCS",mydata$Unique_site_year)
mydata$Unique_site_year=gsub("LHT","LHAT",mydata$Unique_site_year)

# Import corrected dates (I corrected the original database for the FP but not the WL dataset - not sure why)
date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 

newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv")) # note that in the original vegetation database the dates of collection were incorrect for a number of years
newdates$Date.of.collection <- as.Date(newdates$Date.of.collection, format="%d/%m/%Y") 
newdates$Unique_site_year=gsub("KRT","KT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CCNT","NCT",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("CHT","CCS",newdates$Unique_site_year)
newdates$Unique_site_year=gsub("LHT","LHAT",newdates$Unique_site_year)

newdates=newdates[,c(1,3,4,5,6)]
mydata2=merge(mydata, newdates, by="Unique_site_year" ) 

# Assign season to dates

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

