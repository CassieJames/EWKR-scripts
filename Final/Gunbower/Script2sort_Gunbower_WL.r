# Script to turn data into a site by species matrix, 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 2rd June 2018
###########################################################################################
# Import data and check that species have a short species code and its unique

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("Gunbower_WL.csv"))

species=as.data.frame(unique(mydata$Species))
colnames(species)="species.name"

species.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended
mycodes.check=merge(mycodes,species,by.x="Scientific.name",by.y="species.name",all.x=TRUE)

#write.csv(mycodes.check , file = "Gunbower_species_check.csv") # save data out to see which Gunbower codes have not been matched with the master list

mydata=merge(mydata,mycodes,by.x="Species",by.y="Scientific.name",all.x=TRUE)

specieslist=unique(mydata$sp_code_simple)
sitelist=unique(mydata$Icon.Site_Transect_Year_Season)
takeout=c("Inundate", "Lea.litt", "Bar.grou") # remove non species

Output= matrix(NA,nrow=length(sitelist), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=specieslist
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata[which(mydata$Icon.Site_Transect_Year_Season==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code_simple)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]
sitesp=sitesp[!is.na(sitesp)]

if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=max(tdata$Cover..m2.[which(tdata$sp_code_simple==spp)])
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 

rownames(Output)=gsub("Gunbower","GUN",rownames(Output))
rownames(Output)=gsub("Spring","SP",rownames(Output))
rownames(Output)=gsub("Summer","SU",rownames(Output))
rownames(Output)=gsub("Autumn","AU",rownames(Output))
rownames(Output)=gsub("Winter","WI",rownames(Output))

write.csv(Output , file = "Spp_sitetransect_year_season matrix Gunbower_WL_July 2018.csv") # save data out

###################################################################################################################################
#### Summarise by wetland 

GUNsites=c("LL","GS","LG","IP","RL","BLS","FB","CS","LR","COS")

tada = matrix(NA,nrow=length(GUNsites),ncol=ncol(Output))#define the output matrix
rownames(tada)=GUNsites
colnames(tada)=colnames(Output)


for(s in GUNsites) { # 
tdata=Output[grep(s,rownames(Output)),]
out=colSums(tdata)
tada[grep(s,rownames(tada)),] <- out
}

rownames(tada)=paste("GUN_",rownames(tada),sep="")

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_matrix Gunbower_WL_July 2018.csv") # save data out

###################################################################################################################################
#### Summarise by year and wetland 

GUNsites=c("LL","GS","LG","IP","RL","BLS","FB","CS","LR","COS")

years=c(2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2015,2016)

tt = expand.grid(GUNsites,years); tt = paste(tt[,1],"_",tt[,2],sep="")

tada = matrix(NA,nrow=length(tt),ncol=ncol(Output))#define the output matrix
rownames(tada)=tt
colnames(tada)=colnames(Output)
		
for(s in GUNsites) { # 

tdata=Output[grep(s,rownames(Output)),]

for (yy in years) {
ttdata=tdata[grep(yy,rownames(tdata)),]

tada[grep(paste(s,"_",yy,sep=""), rownames(tada), fixed=TRUE),] <-colSums(ttdata)
}}

rownames(tada)=paste("GUN_",rownames(tada),sep="")

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_year matrix Gunbower_WL_July 2018.csv") # save data out

################################################################################################################################
#### Summarise by year and wetland and season

GUNsites=c("LL","GS","LG","IP","RL","BLS","FB","CS","LR","COS")

years=c(2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2015,2016)

seasons=c("_SU", "_SP", "_AU","_WI")

tt = expand.grid(GUNsites,years,seasons); tt = paste("GUN_",tt[,1],"_",tt[,2],tt[,3],sep="")

tada = matrix(NA,nrow=length(tt),ncol=ncol(Output))#define the output matrix
rownames(tada)=tt
colnames(tada)=colnames(Output)
		
for(s in GUNsites) { # 

tdata=Output[grep(s,rownames(Output)),]

for (yy in years) {
ttdata=tdata[grep(yy,rownames(tdata)),]

for (seas in seasons){
tttdata=ttdata[grep(seas,rownames(ttdata)),]

tada[grep(paste(s,"_",yy,seas,sep=""), rownames(tada), fixed=TRUE),] <-colSums(tttdata)
}}}

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records 
tada=tada[,colSums(tada!= 0) > 0]

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_year_season matrix Gunbower_WL_July 2018.csv") # save data out

################################################################################################################################
# Data summary

wetlist=c("LL","GS","LG","IP","RL","BLS","FB","CS","LR","COS")
yoi=c("2000","2001","2002","2003","2004","2005","2006","2007","2008", "2009","2010","2011", "2012", "2013", "2014", "2015","2016")
season=c("SU", "AU","WI","SP")

tt = expand.grid(wetlist); tt = paste("GUN_",tt[,1],sep="")

tada = matrix(NA,nrow=length(tt),ncol=23)#define the output matrix
rownames(tada)=tt
colnames(tada)=c("Totaln","Transectn","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","SU","AU","WI","SP")

for(w in wetlist) { # 

tdata=Output[grep(w,rownames(Output)),]
Totaln=nrow(tdata)
UniTs =length(unique(str_extract(rownames(tdata), "\\d")))
tada[grep(w,rownames(tada)),grep("Totaln",colnames(tada))] = Totaln
tada[grep(w,rownames(tada)),grep("Transectn",colnames(tada))] = UniTs


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

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Data summary table Gunbower_WL_July 2018.csv") # save data out