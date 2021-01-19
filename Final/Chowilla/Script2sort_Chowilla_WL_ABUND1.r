################################################################################################################
#### Chowilla vegetation sorting
#### June 2018
###############################################################################################################
#### Various stages as Chowilla dataset is in two different access databases
#### Step 1 pull Load species names from both database, match with their respective codes, load my master list and match species to my list and extract codes, create single species list
####
###############################################################################################################
#### create single species list 

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydata=data.frame(read.csv("Chowilla_pre2013.csv")) # import pre 2013 data - data is in code form and codes differ from the master list ofcourse !! So need to match Chowilla codes to species and then revert to master code list
chowilla_codes=data.frame(read.csv("Chowilla_Codes_pre2013.csv")) # read in codes
mydata=merge(mydata,chowilla_codes,by="Species",all.x=TRUE) 

species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include chowilla

Chowspecies=as.data.frame(unique(mydata$Taxa))
colnames(Chowspecies)="Taxa"

#mycodes.check=merge(Chowspecies,mycodes,by.x="Taxa",by.y="Scientific.name",all.x=TRUE,all.y=TRUE) # merge data with masterlist code in order to match with codes ion other datasets
#write.csv(mycodes.check , file = "Chowillapre2013_species_check.csv") # save data out to see which  have not been matched with the master list

species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # load amended list 

mydata=merge(mydata,mycodes,by.x="Taxa",by.y="Scientific.name",all.x=TRUE) # merge data with masterlist code in order to match with codes ion other datasets
write.csv(mydata , file = "Chowillapre2013_species_check.csv") # save data out again as some codes don't have species attached

ChowillaTaxaPre2013 <-unique(mydata$sp_code_simple)

###################################################

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydata=data.frame(read.csv("Chowilla_2006-2016.csv")) # import pre 2013 data - data is in code form and codes differ from the master list ofcourse !! So need to match Chowilla codes to species and then revert to master code list
chowilla_codes=data.frame(read.csv("Chowilla_Codes_post2013.csv")) # read in codes
mydata=merge(mydata,chowilla_codes,by.x="Species.Code",by.y="Species",all.x=TRUE) 

species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include chowilla

Chowspecies=as.data.frame(unique(mydata$taxa))
colnames(Chowspecies)="Taxa"

#mycodes.check=merge(Chowspecies,mycodes,by.x="Taxa",by.y="Scientific.name",all.x=TRUE,all.y=TRUE) # merge data with masterlist code in order to match with codes ion other datasets
#write.csv(mycodes.check , file = "Chowillapost2013_species_check.csv") # save data out to see which  have not been matched with the master list

mydata=merge(mydata,mycodes,by.x="taxa",by.y="Scientific.name",all.x=TRUE) # merge data with masterlist code in order to match with codes ion other datasets
ChowillaTaxaPost2013 <-unique(mydata$sp_code_simple)

ChowillaTaxa <-unique(as.character(ChowillaTaxaPre2013),as.character(ChowillaTaxaPost2013)) # create a single species list with correct simplified codes

takeout=c("Inundate", "Lea.litt", "Bar.grou") # remove non species
ChowillaTaxa=ChowillaTaxa[!ChowillaTaxa %in% takeout]

###############################################################################################################
#### Sorting out pre 2013 data into wp_year_season dataset

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydata=data.frame(read.csv("Chowilla_pre2013.csv")) # import pre 2013 data - data is in code form and codes differ from the master list ofcourse !! So need to match Chowilla codes to species and then revert to master code list
chowilla_codes=data.frame(read.csv("Chowilla_Codes_pre2013.csv")) # read in chowilla codes
mydata=merge(mydata,chowilla_codes,by="Species",all.x=TRUE)  # merge data with codes to get full names
species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include chowilla
mydata=merge(mydata,mycodes,by.x="Taxa",by.y="Scientific.name",all.x=TRUE) # merge data with masterlist code in order to match with codes ion other datasets

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))

mydata$Site=as.character(mydata$Site)

wplist=trimws(as.character(unique(mydataGPS$Site.Number)))
wplist=wplist[wplist!=""]# removes empty cells
yoi=c("2004","2005","2006","2007","2008","2009","2010","2011", "2012", "2013", "2014","2015","2016")
seasons=c("Summer", "Spring", "Autumn","Winter")


tt = expand.grid(wplist,yoi, seasons); tt = paste(tt[,1],"_",tt[,2],"_",tt[,3],sep="")

tada = matrix(NA,nrow=length(tt),ncol=length(ChowillaTaxa))#define the output matrix
rownames(tada)=tt
colnames(tada)=ChowillaTaxa

for(w in wplist) { # 
tdata=mydata[which(mydata$Site == w),]

for (yy in yoi) {
ttdata=tdata[grep(yy,tdata$Year),]

for (seas in seasons){

tttdata=ttdata[which(ttdata$Season == seas),]
mysp=unique(tttdata$sp_code_simple)

for (sp in mysp){
tttdata=ttdata[which(ttdata$sp_code_simple==sp),]

wetyr=paste(w,"_",yy,"_",seas,sep="")
tada[grep(wetyr, rownames(tada), fixed=TRUE),grep(sp,colnames(tada), fixed=TRUE)] <-sum(tttdata$Abundance)
}}}}

tada[is.na(tada)]<-0

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)

Site_year_pre2013 <-sort(rownames(tada))

rownames(tada)=gsub("Spring","SP",rownames(tada))
rownames(tada)=gsub("Summer","SU",rownames(tada))
rownames(tada)=gsub("Autumn","AU",rownames(tada))
rownames(tada)=gsub("Winter","WI",rownames(tada))

tada[tada > 0] <- 1

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Spp_wp_year_season matrix Chow2013_May 2019 ABUND.csv") # save data out

###############################################################################################################
#### Sorting out post 2013 data into wp_year_season dataset

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydata=data.frame(read.csv("Chowilla_2006-2016.csv")) # import pre 2013 data - data is in code form and codes differ from the master list ofcourse !! So need to match Chowilla codes to species and then revert to master code list
chowilla_codes=data.frame(read.csv("Chowilla_Codes_post2013.csv")) # read in codes

dates=data.frame(read.csv("Survey dates.csv")) # import dates

left = function(text, num_char) {
  substr(text, 1, num_char)
}

mydata$Date.ID = paste("FP",left(mydata$Survey.ID,6),sep="")
mydata=merge(mydata,dates,by.x="Date.ID",by.y="Site_year",all.x=TRUE,all.y=TRUE) 
mydata$Date <- as.Date(mydata$Date, format="%d/%m/%Y")


library(zoo)
      yq <- as.yearqtr(as.yearmon(mydata$Date , "%Y-%m-%d") + 1/12)
      mydata$Season <- factor(format(yq, "%q"), levels = 1:4, 
      labels = c("SU", "AU", "WI", "SP"))


mydata=merge(mydata,chowilla_codes,by.x="Species.Code",,by.y="Species",all.x=TRUE) 
species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include chowilla
mydata=merge(mydata,mycodes,by.x="taxa",by.y="Scientific.name",all.x=TRUE) # merge data with masterlist code in order to match with codes ion other datasets


data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))

left = function(text, num_char) {
  substr(text, 1, num_char)
}

mid = function(text, start_num, num_char) {
  substr(text, start_num, start_num + num_char - 1)
}

right = function(text, num_char) {
  substr(text, nchar(text) - (num_char-1), nchar(text))
}

mydata$Site <-paste("FP",left(mydata$Survey.ID,3),sep="") # had to create a separate Site column as it was merged with the survey code (site_year_rep)
mydata$Year <-paste("20",mid(mydata$Survey.ID,5,2),sep="")
#mydata$Rep <-paste("rep",right(mydata$Survey.ID,2),sep="") # this does not work so I just added the column to original csv file for ease

wplist=trimws(as.character(unique(mydataGPS$Site.Number)))
wplist=wplist[wplist!=""]# removes empty cells
yoi=c("2004","2005","2006","2007","2008","2009","2010","2011", "2012", "2013", "2014","2015","2016")
seasons=c("SU", "SP", "AU","WI")

tt = expand.grid(wplist,yoi,seasons); tt = paste(tt[,1],"_",tt[,2],"_",tt[,3],sep="")

tada = matrix(NA,nrow=length(tt),ncol=length(ChowillaTaxa))#define the output matrix
rownames(tada)=tt
colnames(tada)=ChowillaTaxa

for(w in wplist) { # 


tdata=mydata[which(mydata$Site == w),]

for (yy in yoi) {
ttdata=tdata[grep(yy,tdata$Year),]

for (seas in seasons){

tttdata=ttdata[which(ttdata$Season == seas),]
mysp=unique(tttdata$sp_code_simple)

for (sp in mysp){
ttttdata=tttdata[which(tttdata$sp_code_simple==sp),]


wetyr=paste(w,"_",yy,"_",seas,sep="")

tada[grep(wetyr, rownames(tada), fixed=TRUE),grep(sp,colnames(tada), fixed=TRUE)] <-sum(ttttdata$Abundance)
}}}}


tada[is.na(tada)]<-0

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)

tada[tada > 0] <- 1

Site_year_post2013 <-sort(rownames(tada))

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Spp_wp_year_season matrix Chow2016_May 2019 ABUND.csv") # save data out

########################################################################################################################################
ChowillaWetlandSites<-c(as.character(Site_year_post2013),as.character(Site_year_pre2013)) # create a single site list by year
Chow2013data=data.frame(read.csv("Spp_wp_year_season matrix Chow2013_May 2019 ABUND.csv",row.names = 1 )) # this list has now been amended to include chowilla
Chow2016data=data.frame(read.csv("Spp_wp_year_season matrix Chow2016_May 2019 ABUND.csv",row.names = 1 )) # this list has now been amended to include chowilla

CHOWALL =rbind(Chow2013data,Chow2016data)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))
mysites=as.character(mydataGPS$Codes)
mysites=unique(mysites[mysites!=""])# removes empty cells
remove=c("MON", "CIL", "BBW")
mysites=mysites [! mysites %in% remove]

png(paste(image.dir,"Chowilla Wetlands accumulation_curves May 2019.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(4,3),cex=1,oma=c(2,0,1,0.5))
		
for(s in mysites) { # 

FPS=mydataGPS[mydataGPS$Codes==s,"Site.Number"]
FPS=trimws(FPS) # removes white spaces

tdata=NULL

for (i in 1:length(FPS)){

tdata=CHOWALL[grep(FPS[i],rownames(CHOWALL)),]

if (i==1) {Tdata=tdata}else {Tdata=rbind(Tdata,tdata)}

}


Tdata[Tdata>0] <-1

Tdata=Tdata[,colSums(Tdata!= 0) > 0]	
Tdata=Tdata[rowSums(Tdata!= 0) > 0,]	

set.seed(101)
accum <- specaccum(Tdata, method="exact")
slopes <- with(accum,diff(richness)/diff(sites))
plat <-which(slopes<0.1)[1]
plateau<-round(accum$richness[plat],0)
Exact<-accum$richness[length(accum$richness)]

sp1 <- specaccum(Tdata, method = "collector")
MM <- fitspecaccum(accum,  "michaelis-menten")
assym <- fitspecaccum(accum,  "asymp")
chao <-round(specpool(Tdata)$chao,0)
boot <-round(specpool(Tdata)$boot,0)

z <- betadiver(Tdata, "z")

plot(accum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlab=s, ylab="Cumulative species number" )
plot(sp1, add=TRUE,col="red")

legend('bottomright', legend = c(paste("Chao =", chao),paste("Exact =",Exact)))

}
dev.off()

#####################################################################################################################################################
##### Chowilla data all together

ChowillaWetlandSites<-c(as.character(Site_year_post2013),as.character(Site_year_pre2013)) # create a single site list by year
Chow2013data=data.frame(read.csv("Spp_wp_year_season matrix Chow2013_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla
Chow2016data=data.frame(read.csv("Spp_wp_year_season matrix Chow2016_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla

CHOWALL =rbind(Chow2013data,Chow2016data)

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))
mysites=as.character(mydataGPS$Codes)
mysites=unique(mysites[mysites!=""])# removes empty cells

tada = matrix(NA,nrow=length(mysites),ncol=length(ChowillaTaxa))#define the output matrix
rownames(tada)=mysites
colnames(tada)=ChowillaTaxa

		
for(s in mysites) { # 

FPS=mydataGPS[mydataGPS$Codes==s,"Site.Number"]
FPS=trimws(FPS) # removes white spaces

tdata=NULL

for (i in 1:length(FPS)){
tdata=CHOWALL[grep(FPS[i],rownames(CHOWALL)),]
if (i==1) {Tdata=tdata}else {Tdata=rbind(Tdata,tdata)}
}

tada[grep(s, rownames(tada), fixed=TRUE),] <-colSums(Tdata)
}

rownames(tada)=paste("CHOW_",rownames(tada),sep="")

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_matrix Chowilla_WL_May 2019.csv") # save data out

#####################################################################################################################################################
##### Chowilla data all together - wetland by year

ChowillaWetlandSites<-c(as.character(Site_year_post2013),as.character(Site_year_pre2013)) # create a single site list by year
Chow2013data=data.frame(read.csv("Spp_wp_year_season matrix Chow2013_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla
Chow2016data=data.frame(read.csv("Spp_wp_year_season matrix Chow2016_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla

CHOWALL =rbind(Chow2013data,Chow2016data)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))
mysites=as.character(mydataGPS$Codes)
mysites=unique(mysites[mysites!=""])# removes empty cells

years=c(2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2015,2016)

tt = expand.grid(mysites,years); tt = paste(tt[,1],"_",tt[,2],sep="")


tada = matrix(NA,nrow=length(mysites)*length(years),ncol=length(ChowillaTaxa))#define the output matrix
rownames(tada)=tt
colnames(tada)=ChowillaTaxa

		
for(s in mysites) { # 

FPS=mydataGPS[mydataGPS$Codes==s,"Site.Number"] # identifies locations within wetlands
FPS=trimws(FPS) # removes white spaces
tdata=CHOWALL[grep(paste(FPS,collapse="|"),rownames(CHOWALL)),]

for (yy in years) {
ttdata=tdata[grep(yy,rownames(tdata)),]

tada[grep(paste(s,"_",yy,sep=""), rownames(tada), fixed=TRUE),] <-colSums(ttdata)
}}

rownames(tada)=paste("CHOW_",rownames(tada),sep="")

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_year matrix Chowilla_WL_May 2019.csv") # save data out

#####################################################################################################################################################
##### Chowilla data all together - wetland by year by season

ChowillaWetlandSites<-c(as.character(Site_year_post2013),as.character(Site_year_pre2013)) # create a single site list by year
Chow2013data=data.frame(read.csv("Spp_wp_year_season matrix Chow2013_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla
Chow2016data=data.frame(read.csv("Spp_wp_year_season matrix Chow2016_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla

CHOWALL =rbind(Chow2013data,Chow2016data)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))
mysites=as.character(mydataGPS$Codes)
mysites=unique(mysites[mysites!=""])# removes empty cells

years=c(2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2015,2016)
seasons=c("_SU", "_SP", "_AU","_WI")

tt = expand.grid(mysites,years,seasons); tt = paste("CHOW_",tt[,1],"_",tt[,2],tt[,3],sep="")


tada = matrix(NA,nrow=length(mysites)*length(years)*length(seasons),ncol=length(ChowillaTaxa))#define the output matrix
rownames(tada)=tt
colnames(tada)=ChowillaTaxa

		
for(s in mysites) { # 

FPS=mydataGPS[mydataGPS$Codes==s,"Site.Number"] # identifies locations within wetlands
FPS=trimws(FPS) # removes white spaces
tdata=CHOWALL[grep(paste(FPS,collapse="|"),rownames(CHOWALL)),]

for (yy in years) {
ttdata=tdata[grep(yy,rownames(tdata)),]

tada[grep(paste(s,"_",yy,sep=""), rownames(tada), fixed=TRUE),] <-colSums(ttdata)
}}


data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_year_season matrix Chowilla_WL_May 2019.csv") # save data out

###################################################################################################################################
# May 2019 - this creates the dataframe to use for the common species ONLY beta diversity analysis
###################################################################################################################################


Output=ceiling(tada) # first round up as these numbers represent % cover - any detection is a detection
Output[Output > 0] <- 1# Second change to presence/absence because we are just interested in repeated observations over time

CHOWsites=mysites

tada = matrix(NA,nrow=length(CHOWsites),ncol=ncol(Output))#define the output matrix
rownames(tada)=CHOWsites
colnames(tada)=colnames(Output)

for(s in CHOWsites) { # 
tdata=Output[grep(s,rownames(Output)),]
out=colSums(tdata)
tada[grep(s,rownames(tada)),] <- out
}

rownames(tada)=paste("CHOW","_",rownames(tada),sep="")

tada[tada<3 ] <- 0 # if species recorded less than 5 times turn value to zero

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_matrix Chowilla_WL_May 2019 COMMON ONLY.csv") # save data out


#####################################################################################################################################################
##### Chowilla data summary

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))
mysites=as.character(mydataGPS$Codes)
wetlist=unique(mysites[mysites!=""])# removes empty cells

ChowillaWetlandSites<-c(as.character(Site_year_post2013),as.character(Site_year_pre2013)) # create a single site list by year
Chow2013data=data.frame(read.csv("Spp_wp_year_season matrix Chow2013_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla
Chow2016data=data.frame(read.csv("Spp_wp_year_season matrix Chow2016_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla
CHOWALL =rbind(Chow2013data,Chow2016data)

yoi=c("2000","2001","2002","2003","2004","2005","2006","2007","2008", "2009","2010","2011", "2012", "2013", "2014", "2015","2016")
season=c("SU", "AU","WI","SP")

tt = expand.grid(wetlist); tt = paste("CHOW_",tt[,1],sep="")

tada = matrix(NA,nrow=length(tt),ncol=24)#define the output matrix
rownames(tada)=tt
colnames(tada)=c("Totaln","Transectn","Elevn","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","SU","AU","WI","SP")
library(stringr)

for(w in wetlist) { # 

FPS=mydataGPS[mydataGPS$Codes==w,"Site.Number"] # identifies locations within wetlands
FPS=trimws(FPS) # removes white spaces
tdata=CHOWALL[grep(paste(FPS,collapse="|"),rownames(CHOWALL)),]
Totaln=nrow(tdata)
tada[grep(w,rownames(tada)),grep("Totaln",colnames(tada))] = Totaln

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
data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Data summary table Chowilla_WL_July 2018.csv") # save data out


###############################################################################################################
#### Including replicates as separate lines pre 2013 data into wp_rep_year_season dataset

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydata=data.frame(read.csv("Chowilla_pre2013.csv")) # import pre 2013 data - data is in code form and codes differ from the master list ofcourse !! So need to match Chowilla codes to species and then revert to master code list
chowilla_codes=data.frame(read.csv("Chowilla_Codes_pre2013.csv")) # read in chowilla codes
mydata=merge(mydata,chowilla_codes,by="Species",all.x=TRUE)  # merge data with codes to get full names
species.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include chowilla
mydata=merge(mydata,mycodes,by.x="Taxa",by.y="Scientific.name",all.x=TRUE) # merge data with masterlist code in order to match with codes ion other datasets

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))

mydata$Site=as.character(mydata$Site)

wplist=trimws(as.character(unique(mydataGPS$Site.Number)))
wplist=wplist[wplist!=""]# removes empty cells
yoi=c("2004","2005","2006","2007","2008","2009","2010","2011", "2012", "2013", "2014","2015","2016")
seasons=c("Summer", "Spring", "Autumn","Winter")
rep=c(1:3)

tt = expand.grid(wplist,c(1:3),yoi, seasons); tt = paste(tt[,1],"_",tt[,2],"_",tt[,3],"_",tt[,4],sep="")

tada = matrix(NA,nrow=length(tt),ncol=length(ChowillaTaxa))#define the output matrix
rownames(tada)=tt
colnames(tada)=ChowillaTaxa

for(w in wplist) { # 
tdata=mydata[which(mydata$Site == w),]

for(r in rep) { # 
ttdata=tdata[which(tdata$Rep == r),]

for (yy in yoi) {
tttdata=ttdata[grep(yy,ttdata$Year),]

for (seas in seasons){

ttttdata=tttdata[which(tttdata$Season == seas),]
mysp=unique(ttttdata$sp_code_simple)

for (sp in mysp){
spdata=ttttdata[which(ttttdata$sp_code_simple==sp),]
wetyr=paste(w,"_",r,"_",yy,"_",seas,sep="")
tada[grep(wetyr, rownames(tada), fixed=TRUE),grep(sp,colnames(tada), fixed=TRUE)] <-sum(spdata$Abundance)
}}}}}

tada[is.na(tada)]<-0

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)

Site_year_pre2013 <-sort(rownames(tada))

rownames(tada)=gsub("Spring","SP",rownames(tada))
rownames(tada)=gsub("Summer","SU",rownames(tada))
rownames(tada)=gsub("Autumn","AU",rownames(tada))
rownames(tada)=gsub("Winter","WI",rownames(tada))

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Spp_wp_rep_year_season matrix Chow2013_July 2018.csv") # save data out

###############################################################################################################
#### Sorting out post 2013 data into wp_rep_year_season dataset

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydata=data.frame(read.csv("Chowilla_2006-2016.csv")) # import pre 2013 data - data is in code form and codes differ from the master list ofcourse !! So need to match Chowilla codes to species and then revert to master code list
chowilla_codes=data.frame(read.csv("Chowilla_Codes_post2013.csv")) # read in codes

dates=data.frame(read.csv("Survey dates.csv")) # import dates

left = function(text, num_char) {
  substr(text, 1, num_char)
}

mydata$Date.ID = paste("FP",left(mydata$Survey.ID,6),sep="")
mydata=merge(mydata,dates,by.x="Date.ID",by.y="Site_year",all.x=TRUE,all.y=TRUE) 
mydata$Date <- as.Date(mydata$Date, format="%d/%m/%Y")

#new_mydata <- mydata[rowSums(is.na(mydata)) > 0,] # returns all rows with NA in them to check which dates are missing - have emailed Jason
#write.csv(new_mydata , file = "Chowillapost2013 with date issues.csv") # save data out again as some codes don't have species attached

library(zoo)
      yq <- as.yearqtr(as.yearmon(mydata$Date , "%Y-%m-%d") + 1/12)
      mydata$Season <- factor(format(yq, "%q"), levels = 1:4, 
      labels = c("SU", "AU", "WI", "SP"))


mydata=merge(mydata,chowilla_codes,by.x="Species.Code",,by.y="Species",all.x=TRUE) 
species.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to include chowilla
mydata=merge(mydata,mycodes,by.x="taxa",by.y="Scientific.name",all.x=TRUE) # merge data with masterlist code in order to match with codes ion other datasets


data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))

left = function(text, num_char) {
  substr(text, 1, num_char)
}

mid = function(text, start_num, num_char) {
  substr(text, start_num, start_num + num_char - 1)
}

right = function(text, num_char) {
  substr(text, nchar(text) - (num_char-1), nchar(text))
}

mydata$Site <-paste("FP",left(mydata$Survey.ID,3),sep="") # had to create a separate Site column as it was merged with the survey code (site_year_rep)
mydata$Year <-paste("20",mid(mydata$Survey.ID,5,2),sep="")
#mydata$Rep <-paste("rep",right(mydata$Survey.ID,2),sep="") # this does not work so I just added the column to original csv file for ease

wplist=trimws(as.character(unique(mydataGPS$Site.Number)))
wplist=wplist[wplist!=""]# removes empty cells
yoi=c("2004","2005","2006","2007","2008","2009","2010","2011", "2012", "2013", "2014","2015","2016")
seasons=c("SU", "SP", "AU","WI")
rep=c(1:3)

tt = expand.grid(wplist,rep,yoi,seasons); tt = paste(tt[,1],"_",tt[,2],"_",tt[,3],"_",tt[,4],sep="")

tada = matrix(NA,nrow=length(tt),ncol=length(ChowillaTaxa))#define the output matrix
rownames(tada)=tt
colnames(tada)=ChowillaTaxa

for(w in wplist) { # 
tdata=mydata[which(mydata$Site == w),]

for(r in rep) { # 
ttdata=tdata[which(tdata$REP == r),]

for (yy in yoi) {
tttdata=ttdata[grep(yy,ttdata$Year),]

for (seas in seasons){

ttttdata=tttdata[which(tttdata$Season == seas),]
mysp=unique(ttttdata$sp_code_simple)

for (sp in mysp){
spdata=ttttdata[which(ttttdata$sp_code_simple==sp),]


wetyr=paste(w,"_",r,"_",yy,"_",seas,sep="")

tada[grep(wetyr, rownames(tada), fixed=TRUE),grep(sp,colnames(tada), fixed=TRUE)] <-sum(spdata$Abundance)
}}}}}


tada[is.na(tada)]<-0

tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records (e.g. BIT was not sampled until 2013 so there is no data for 2008-2012)


Site_year_post2013 <-sort(rownames(tada))

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Spp_wp_rep_year_season matrix Chow2016_July 2018.csv") # save data out

########################################################################################################################################
ChowillaWetlandSites<-c(as.character(Site_year_post2013),as.character(Site_year_pre2013)) # create a single site list by year
Chow2013data=data.frame(read.csv("Spp_wp_rep_year_season matrix Chow2013_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla
Chow2016data=data.frame(read.csv("Spp_wp_rep_year_season matrix Chow2016_July 2018.csv",row.names = 1 )) # this list has now been amended to include chowilla

CHOWALLREPS =rbind(Chow2013data,Chow2016data)

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd (data.dir)
mydataGPS=data.frame(read.csv("Chowilla_GPS.csv"))
mysites=as.character(mydataGPS$Codes)
mysites=unique(mysites[mysites!=""])# removes empty cells
remove=c("MON", "CIL", "BBW")
mysites=mysites [! mysites %in% remove]

png(paste(image.dir,"Chowilla Wetlands accumulation_curves_with replicates sep PA.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(4,3),cex=1,oma=c(2,0,1,0.5))
		
for(s in mysites) { # 

FPS=mydataGPS[mydataGPS$Codes==s,"Site.Number"]
FPS=trimws(FPS) # removes white spaces

tdata=NULL

for (i in 1:length(FPS)){

tdata=CHOWALLREPS[grep(FPS[i],rownames(CHOWALLREPS)),]

if (i==1) {Tdata=tdata}else {Tdata=rbind(Tdata,tdata)}

}


Tdata[Tdata>0] <-1

Tdata=Tdata[,colSums(Tdata!= 0) > 0]	
Tdata=Tdata[rowSums(Tdata!= 0) > 0,]	

set.seed(101)
accum <- specaccum(Tdata, method="exact")
slopes <- with(accum,diff(richness)/diff(sites))
plat <-which(slopes<0.1)[1]
plateau<-round(accum$richness[plat],0)
Exact<-accum$richness[length(accum$richness)]

sp1 <- specaccum(Tdata, method = "collector")
MM <- fitspecaccum(accum,  "michaelis-menten")
assym <- fitspecaccum(accum,  "asymp")
chao <-round(specpool(Tdata)$chao,0)
boot <-round(specpool(Tdata)$boot,0)

z <- betadiver(Tdata, "z")

plot(accum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlab=s, ylab="Cumulative species number" )
plot(sp1, add=TRUE,col="red")

legend('bottomright', legend = c(paste("Chao =", chao),paste("Exact =",Exact)))

}
dev.off()
