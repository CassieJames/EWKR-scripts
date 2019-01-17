##########################################################################################################
### Script to plot modelled water depths against observations to validate water depth (as far as is possible!)
### C James September 2018

library(hydrostats)
library(tidyr)
library(ggplot2)

# step 1 - upload all required data 

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Depth.dat=read.csv("Hattah_sites_time_series_dataV3.csv", check.names=FALSE) # load hydrological information but prevent r form converting column names to a more 'friendly' format
Depth_long <- gather(Depth.dat, site, depth, gather_cols=2:406, factor_key=TRUE) # had to use numbers for columns as function would not recognise names

Depth_long$site <-as.character(Depth_long$site)
temp <- strsplit(Depth_long$site, split="_")
SITE_only <- sapply(temp, "[", 1)
Depth_long$SITE <-SITE_only # this dataset has both ends of the transect - will need to use grep to match with vegetation dataset which does not include a or b
Depth_long$DATE=as.Date(as.character(Depth_long$DATE), format = "%m/%d/%Y") # R needs to recognise column as dates


# Bring in wetland data
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_WL_FGcorrections.csv"))


# Bring in corrected sample dates and merge with mydata
date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 

newdates=data.frame(read.csv("HTH_WL_dates_corrected.csv")) # note that in the original vegetation database the dates of collection were incorrect for a number of years
newdates$Date.of.collection <- as.Date(newdates$Date.of.collection, format="%d/%m/%Y") 
mydata <- merge(newdates, mydata,by="Unique_site_year")

Sites_sampled <-unique(mydata[,c("Site.ID.x","Date.of.collection.x")]) # creates unique list of sites and dates sampled - with corrected 
Sites_sampled$Date.of.collection<-as.Date(as.character(Sites_sampled$Date.of.collection.x), format = "%d/%m/%Y")
Sites_sampled$Site.ID.x <- as.character(Sites_sampled$Site.ID.x)
colnames(Sites_sampled)=c("Site.ID", "Date.of.collection")

Hydro_Sites_Sampled=Sites_sampled
Hydro_Sites_Sampled$DOD <-0 # Depth On Day of sampling


for (i in 1:nrow(Sites_sampled)){ # 
soi=Sites_sampled[i,1]
soi=ifelse(grepl("CCNT",soi),gsub("CCNT","NCT",soi),soi) # this adjusts for the differences in labels between the veg database and the location file sent by Cherie
soi=ifelse(grepl("KRT",soi),gsub("KRT","KT",soi),soi) 
DepOI = as.data.frame(Depth.dat[,grepl(soi, colnames(Depth.dat), fixed=TRUE)])
DepOI=cbind(Depth.dat$DATE,DepOI)
colnames(DepOI)[1]="Date"
DepOI$Date=as.Date(as.character(DepOI$Date), format = "%m/%d/%Y") 
doi=Sites_sampled[i,2]
numcols=ncol(DepOI)
Depth_on_day <-DepOI[which(DepOI$Date==doi),2:numcols]
Hydro_Sites_Sampled[i,c("DOD")] <-ifelse(numcols==2,Depth_on_day,rowMeans(Depth_on_day)) 
}

tada<-Hydro_Sites_Sampled[order(as.Date(Hydro_Sites_Sampled$Date.of.collection, format="%d/%m/%Y")),]
Depth.dat$DATE=as.Date(as.character(Depth.dat$DATE), format = "%m/%d/%Y")

#########################################################################################################################
#### Plots to look at long term depth records for wetlands and plots 
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
Depth.dat=read.csv("Hattah_sites_time_series_dataV3.csv", check.names=FALSE) # load hydrological information but prevent r form converting column
Depth.dat$DATE=as.Date(as.character(Depth.dat$DATE), format = "%d/%m/%Y")

plot.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/WL Hydrology plots/"
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) # set working directory
obs.dat=read.csv("Inundation_obs.csv", check.names=FALSE) 
colnames(obs.dat)[1]="site"

Sites=c("BIT1","BIT2","BIT3","BIT4", "BOT1","BOT2","BOT3","BRT1","BRT2","BRT3","BLT1","BLT2","BLT3","BLT4","CHT1+","CHT2+","CHT3+","CHT4+","CHT1S", 
"CHT3S","CHT2N","CHT4N","NCT1","NCT2","NCT3","NCT4", "HT1","HT2","HT3","HT4", "KT1","KT2","KT3","KT4", "LHT1","LHT2","LHT3",
 "LHT4","MOT1","MOT2","MOT3","MOT4", "NN1","NN2","NN3","NN4", "YT1","YT2","YT3","YT4")

for (s in Sites) {

Hattah.dat = as.data.frame(Depth.dat[,grepl(paste("^",s,sep=""), colnames(Depth.dat))])
image_size=(ceiling((ncol(Hattah.dat))/3))*6.25 # adjusts height of image depending upon number of images

obs.dat$site <-gsub("CCNT", "NCT", obs.dat$site)
obs.dat$site <-gsub("KRT", "KT", obs.dat$site)

Hattah.obs.dat =obs.dat[grep(paste("^",s,sep=""),obs.dat$site),]
Inundated=Hattah.obs.dat[which(Hattah.obs.dat$Inundated=="TRUE"),]
NotInundated=Hattah.obs.dat[which(Hattah.obs.dat$Inundated=="FALSE"),]

Inundated$Date=as.Date(as.character(Inundated$Date), format = "%d/%m/%Y") # dates when field obs suggest site was inundated
NotInundated$Date=as.Date(as.character(NotInundated$Date), format = "%d/%m/%Y") # dates when field obs suggest site was inundated

Hattah.dat <- cbind(Depth.dat$DATE, Hattah.dat)
colnames(Hattah.dat)[1] <-"Date"
Hattah_long <- gather(Hattah.dat, site, depth, -Date, factor_key=TRUE)
Hattah_long$site <-gsub("_.*", "", Hattah_long$site)
Hattah_long$site<-sub("B$", '', Hattah_long$site)
Hattah_long$site<-sub("A$", '', Hattah_long$site)
Hattah_long$site<-as.factor(Hattah_long$site)
Hattah_long$site<-as.factor(Hattah_long$site)
Hattah_agg <- aggregate(Hattah_long$depth,  by = list(Hattah_long$site, Hattah_long$Date), mean)
colnames(Hattah_agg) <-c("site", "Date", "Depth")

png(paste(plot.dir,"/",s,"v5.png",sep=''),width=25, height=image_size, units='cm', res=300, pointsize=15, bg='white')


plt<-ggplot(Hattah_agg, aes(x = Date, y = Depth)) + geom_line(color = "dodgerblue3",fill="dodgerblue3") + ylab("Height (m) ")+facet_wrap(~site,scales = "free_x", ncol=3)+scale_x_date(date_breaks = "1 year")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
plt<-plt+geom_point(mapping=aes(x=Date, y=-0.25),colour='red',data=NotInundated)+facet_wrap(~site,scales = "free_x", ncol=3)
plt<-plt+geom_point(mapping=aes(x=Date, y=-0.25),colour='blue',data=Inundated)+facet_wrap(~site,scales = "free_x", ncol=3)
print(plt)
dev.off()

}