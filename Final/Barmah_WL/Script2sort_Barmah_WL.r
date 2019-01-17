# Script to turn data into a site by species matrix, 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 2rd June 2018
###########################################################################################
# Import data and check that species have a short species code and its unique

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("Copy of BMFVEGCLEANED.csv"))

#species=as.data.frame(unique(mydata$species_))
#colnames(species)="species.name"

species.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to all PFG added

#mycodes.check=merge(mycodes,species,by.x="Scientific.name",by.y="species.name",all.x=TRUE,all.y=TRUE)
#write.csv(mycodes.check , file = "Barmah_species_check.csv") # save data out to see which codes have not been matched with the master list

mydata=merge(mydata,mycodes,by.x="species_",by.y="Scientific.name",all.x=TRUE)

specieslist=unique(mydata$sp_code_simple) # take species codes to create column headings for site by species matirx
specieslist=specieslist[!is.na(specieslist)]

sitelist=unique(mydata$Unique_id)

Output= matrix(NA,nrow=length(sitelist), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=specieslist
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata[which(mydata$Unique_id==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code_simple)
sitesp=sitesp[!is.na(sitesp)]

if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=sum(tdata$cover_[which(tdata$sp_code_simple==spp & tdata$Vartype =="A")]) # I have summed across transects for the time being...
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix (matrix is BARM_SITE_YEAR_SEASON)

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir) 

write.csv(Output , file = "Spp_site_year_season matrix Barmah_WL_July2018.csv") # save data out

###################################################################################################################################
#### Summarise by wetland 

Barmsites=c("WL","BDEAD","SP","TL","DUCK","RBS","TIB","BG","Alga","TIO","LRS")

tada = matrix(NA,nrow=length(Barmsites),ncol=ncol(Output))#define the output matrix
rownames(tada)=Barmsites
colnames(tada)=colnames(Output)

for(s in Barmsites) { # 
tdata=Output[grep(s,rownames(Output)),]
out=colSums(tdata)
tada[grep(s,rownames(tada)),] <- out
}

rownames(tada)=paste("BARM","_",rownames(tada),sep="")

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_matrix Barmah_WL_July 2018.csv") # save data out

###################################################################################################################################
#### Summarise by year and wetland 

Barmsites=c("WL","BDEAD","SP","TL","DUCK","RBS","TIB","BG","Alga","TIO","LRS")

years=c(1990,1991,1992,1993,1994,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017) #o

tt = expand.grid(Barmsites,years); tt = paste(tt[,1],"_",tt[,2],sep="")

tada = matrix(NA,nrow=length(tt),ncol=ncol(Output))#define the output matrix
rownames(tada)=tt
colnames(tada)=colnames(Output)
		
for(s in Barmsites) { # 

tdata=Output[grep(s,rownames(Output)),]

for (yy in years) {
ttdata=tdata[grep(yy,rownames(tdata)),]

tada[grep(paste(s,"_",yy,sep=""), rownames(tada), fixed=TRUE),] <-colSums(ttdata)
}}

rownames(tada)=paste("BARM","_",rownames(tada),sep="")

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir) 
write.csv(tada , file = "Spp_site_year matrix Barmah_WL_July 2018.csv") # save data out

################################################################################################################################
#### Summarise by wetland, replicate, year and season

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("Copy of BMFVEGCLEANED.csv"))

species.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended to all PFG added

mydata=merge(mydata,mycodes,by.x="species_",by.y="Scientific.name",all.x=TRUE)

specieslist=unique(mydata$sp_code_simple) # take species codes to create column headings for site by species matirx
specieslist=specieslist[!is.na(specieslist)]

sitelist=unique(mydata$Unique_id_trans)

Output= matrix(NA,nrow=length(sitelist), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=specieslist
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata[which(mydata$Unique_id_trans==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code_simple)
sitesp=sitesp[!is.na(sitesp)]

if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=sum(tdata$cover_[which(tdata$sp_code_simple==spp & tdata$Vartype =="A")]) # I have summed across transects for the time being...
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix 

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir) 

write.csv(Output , file = "Spp_site_trans_year_season matrix Barmah_WL_July2018.csv") # save data out

################################################################################################################################
# Data summary

wetlist=c("WL","BDEAD","SP_","TL","DUCK","RBS","TIB","BG","Alga","TIO","LRS")
yoi=c("1990","1991","1992", "1993", "1994", "1995","1996","2000","2001","2002","2003","2004","2005","2006","2007","2008", "2009","2010","2011", "2012", "2013", "2014", "2015","2016")
season=c("SU", "AU","WI","SP")
transect=c("_1_", "_2_", "_3_")

tt = expand.grid(wetlist); tt = paste("Barm_",tt[,1],sep="")

tada = matrix(NA,nrow=length(tt),ncol=30)#define the output matrix
rownames(tada)=tt
colnames(tada)=c("Totaln","Transectn","1990","1991","1992", "1993", "1994", "1995","1996","2000","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","SU","AU","WI","SP")

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

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir)
write.csv(tada , file = "Data summary table Barmah_WL_July 2018.csv") # save data out

################################################################################################################################
# Barmah species accumulation curves
image.dir = "C:/Users/jc246980/Documents/MD Vegetation/Plots/"
wetlist=c("WL","BDEAD","SP_","TL","DUCK","RBS","TIB","BG","Alga","TIO","LRS")

png(paste(image.dir,"Barmah Wetlands accumulation_curves_with replicates sep PA.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(4,3),cex=1,oma=c(2,0,1,0.5))
		
for(w in wetlist) { # 

Tdata=Output[grep(w,rownames(Output)),]
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



plot(accum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlab=w, ylab="Cumulative species number" )
plot(sp1, add=TRUE,col="red")

legend('bottomright', legend = c(paste("Chao =", chao),paste("Exact =",Exact)))

}
dev.off()


