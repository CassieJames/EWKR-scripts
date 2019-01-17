###########################################################################################
# Script to turn data into a site by species matrix
# Part 1 summarises data by wetland and year
# Part 2 summarises by transect and year
# GNU General Public License .. feel free to use / distribute ... no warranties
# 15th December 2017
###########################################################################################
# LMW data summarised to wetland by year

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("LMW_WL.csv"))
species=unique(mydata$Scientific.name)

species.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include all sites
mydata=merge(mydata,mycodes,by="Scientific.name")

#write.csv(mydata , file = "Species code check matrix LMW_WL_August 2018.csv") # save data out to check list

specieslist=unique(mydata$sp_code_simple)
sitelist=unique(mydata$Site_year)
takeout=c("Inundate", "Lea.litt", "Bar.grou") # remove non species
specieslist=specieslist[!specieslist %in% takeout]
Output= matrix(NA,nrow=length(unique(mydata$Site_year)), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=sort(specieslist)
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata[which(mydata$Site_year==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code_simple)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]
if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$sp_code_simple==spp)])
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix

rownames(Output)=paste("LMW_",rownames(Output),sep="")
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


data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 
write.csv(Output , file = "Spp_site_year matrix LMW_WL_July 2018.csv") # save data out

###########################################################################################
# Summarise LMW to site

wetlist=c("BB","CR","LP","UL","MUH","BI","UMWC","W33","SCB","MLH","WL","WW")

tada = matrix(NA,nrow=length(wetlist),ncol=ncol(Output))#define the output matrix
rownames(tada)=wetlist
colnames(tada)=colnames(Output)

for(w in wetlist) { # 
tdata=Output [grep(w,rownames(Output)),]
tdata=colSums(tdata)

tada[grep(w,rownames(tada)),]<- tdata

}

tada[is.na(tada)]<-0

rownames(tada)=paste("LMW_",rownames(tada),sep="")

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site matrix LMW_WL_July 2018.csv") # save data out

###########################################################################################
# LMW data summarised to wetland by year and season

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("LMW_WL.csv"))
species=unique(mydata$Scientific.name)

species.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAH.csv")) # this list has now been amended to include all sites
mydata=merge(mydata,mycodes,by="Scientific.name")

mydata$Date.of.collection <- as.Date(mydata$Date.of.collection, format="%m/%d/%Y")

library(zoo)
      yq <- as.yearqtr(as.yearmon(mydata$Date.of.collection , "%m/%d/%Y") + 1/12)
      mydata$Season <- factor(format(yq, "%q"), levels = 1:4, 
      labels = c("SU", "AU", "WI", "SP"))


mydata$Unique_site_year_season = paste(mydata$Site_year,"_",mydata$Season,sep="")
 
specieslist=unique(mydata$sp_code_simple)
sitelist=unique(mydata$Unique_site_year_season)
takeout=c("Inundate", "Lea.litt", "Bar.grou") # remove non species
specieslist=specieslist[!specieslist %in% takeout]
Output= matrix(NA,nrow=length(unique(mydata$Unique_site_year_season)), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=sort(specieslist)
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata[which(mydata$Unique_site_year_season==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code_simple)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]
if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$sp_code_simple==spp)])
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix

rownames(Output)=paste("LMW_",rownames(Output),sep="")
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

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 
write.csv(Output , file = "Spp_site_year_season matrix LMW_WL_July 2018.csv") # save data out

#######################################################################################################
# LMW data site by transect and year and season

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir)
mydata=data.frame(read.csv("LMW_WL.csv"))
species.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAH.csv")) # this list has now been amended to include chowilla

mid = function(text, start_num, num_char) {
  substr(text, start_num, start_num + num_char - 1)
}

mydata$Unique_site_id <-paste(mydata$Site.ID,"_",mydata$Short.year,sep="")

mydata=merge(mydata,mycodes,by="Scientific.name")
mydata <- within(mydata, Date.of.collection <- as.Date(as.character(Date.of.collection), format = "%m/%d/%Y")) # ensure dates are recognised as dates!

mydata$Date.of.collection <- as.Date(mydata$Date.of.collection, format="%m/%d/%Y")

library(zoo)
      yq <- as.yearqtr(as.yearmon(mydata$Date.of.collection , "%m/%d/%Y") + 1/12)
      mydata$Season <- factor(format(yq, "%q"), levels = 1:4, 
      labels = c("SU", "AU", "WI", "SP"))

# note that short site names do not match site.id column in original database so having to fix up here
mydata$Unique_site_id=gsub("L2Vof","SCB",mydata$Unique_site_id)
	  
mydata$Unique_site_year_season = paste(mydata$Unique_site_id,"_",mydata$Season,sep="")

# Create species x site matrix
specieslist=unique(mydata$sp_code_simple)
sitelist=unique(mydata$Unique_site_year_season)
takeout=c("Inundate", "Lea.litt", "Bar.grou") # remove non species
specieslist=specieslist[!specieslist %in% takeout]
Output= matrix(NA,nrow=length(unique(mydata$Unique_site_year_season)), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=sort(specieslist)
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata[which(mydata$Unique_site_year_season==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code_simple)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]

if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$sp_code_simple==spp)])
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 

rownames(Output)=paste("LMW_",rownames(Output),sep="")
rownames(Output)=gsub("_13","_2013",rownames(Output))
rownames(Output)=gsub("_16","_2016",rownames(Output))
rownames(Output)=gsub("_15","_2015",rownames(Output))
rownames(Output)=gsub("_14","_2014",rownames(Output))
rownames(Output)=gsub("_12","_2012",rownames(Output))
rownames(Output)=gsub("_11","_2011",rownames(Output))
rownames(Output)=gsub("_10","_2010",rownames(Output))
rownames(Output)=gsub("_9","_2009",rownames(Output))
rownames(Output)=gsub("_8","_2008",rownames(Output))
rownames(Output)=gsub("_7","_2007",rownames(Output))
rownames(Output)=gsub("_6","_2006",rownames(Output))

# note that short site names do not match site.id column in original database so having to fix up here
rownames(Output)=gsub("_MLT","_MLH",rownames(Output)) # this code has to be changed otherwise LT is a subset of MLT
rownames(Output)=gsub("_LT","_LP",rownames(Output))
rownames(Output)=gsub("_MUT","_MUH",rownames(Output))

write.csv(Output , file = "Spp_site_transect_year_transect matrix LMW_July 2018.csv") # save data out

##########################################################################################################################
#### Script to summarise data for each wetland.

wetlist=c("BB","CR","LP","UL","MUH","BI","UMWC","W33","SCB","MLH","WL","WW")
yoi=c("2000","2001","2002","2003","2004","2005","2006","2007","2008", "2009","2010","2011", "2012", "2013", "2014", "2015","2016")
season=c("SU", "AU","WI","SP")

tt = expand.grid(wetlist); tt = paste("LMW_",tt[,1],sep="")

tada = matrix(NA,nrow=length(tt),ncol=24)#define the output matrix
rownames(tada)=tt
colnames(tada)=c("Totaln","Transectn","Elevn","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","SU","AU","WI","SP")
library(stringr)

for(w in wetlist) { # 

tdata=Output[grep(w,rownames(Output)),]
Totaln=nrow(tdata)
Units=vapply(strsplit(rownames(tdata),"\\+"), `[`, 1, FUN.VALUE=character(1))
UniTs =length(unique(str_extract(Units, "\\d$")))

Elevs=vapply(strsplit(rownames(tdata),"\\+"), `[`, 2, FUN.VALUE=character(1)) # find any positive elevations
Elevs2=unique(vapply(strsplit(Elevs,"\\_"), `[`, 1, FUN.VALUE=character(1))) # extract number of positive different elevations for site before "_"
Elevs2=Elevs2[!is.na(Elevs2)] # remove Nas
Elevs2=length(unique(str_extract(Elevs2, "\\d+"))) # extract unique numerical values ignoring N S E or W


Elevsminus=vapply(strsplit(rownames(tdata),"\\-"), `[`, 2, FUN.VALUE=character(1)) # find any minus elevations
Elevsminus2=unique(vapply(strsplit(Elevsminus,"\\_"), `[`, 1, FUN.VALUE=character(1))) # extract number of different elevations for site that are negative
Elevsminus2=Elevsminus2[!is.na(Elevsminus2)]
Elevsminus2=length(unique(str_extract(Elevsminus2, "\\d+"))) # extract unique numerical values ignoring N S E or W
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

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Data summary table LMW_WL_July 2018.csv") # save data out