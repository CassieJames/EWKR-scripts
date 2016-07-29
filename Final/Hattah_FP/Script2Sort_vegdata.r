# Script to turn data into a site by species matrix


data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 

mydata=data.frame(read.csv("HTH_FP.csv"))

specieslist=unique(mydata$Scientific_name)
sitelist=unique(mydata$Unique_site_year)

Output= matrix(NA,nrow=length(unique(mydata$Unique_site_year)), ncol=length(unique(mydata$Scientific_name)))

rownames(Output)=unique(mydata$Unique_site_year)


colnames(Output)=sort(unique(mydata$Scientific_name))

for(s in sitelist) {
tdata=mydata[which(mydata$Unique_site_year==s),]
sitesp=unique(tdata$Scientific_name)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$Scientific_name==spp)])
Output[grep(s, rownames(Output)),grep(spp,colnames(Output))] <- abund
}
}

Output[is.na(Output)] <- 0


# Create species matrix with matched environmental data

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 

mydata=data.frame(read.csv("HTH_FP.csv"))

image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/"; setwd(image.dir)

Rainfall.dat=data.frame(read.csv("Rainfall_HTH_FP.csv"))
Rainfall.dat <- within(Rainfall.dat, Date <- as.Date(as.character(Date), format = "%d/%m/%Y"))

mydata <- within(mydata, Date.of.collection <- as.Date(as.character(Date.of.collection), format = "%d/%m/%Y"))

Rainfall.dat=Rainfall.dat[!duplicated(Rainfall.dat), ]
mydata_env=merge(mydata, Rainfall.dat, by.x="Date.of.collection", by.y="Date")

mydata_rainfall=mydata_env[,c("Unique_site_year", "d30", "d90", "d180", "d365")]
mydata_rainfall=unique(mydata_rainfall)

specieslist=unique(mydata$Scientific_name)
sitelist=unique(mydata$Unique_site_year)

Output= matrix(NA,nrow=length(unique(mydata$Unique_site_year)), ncol=length(unique(mydata$Scientific_name)))

rownames(Output)=unique(mydata$Unique_site_year)
colnames(Output)=sort(unique(mydata$Scientific_name))

Output=as.data.frame(Output)
Output=merge(Output, mydata_rainfall, by.x="row.names", by.y="Unique_site_year", all.x=TRUE)

for(s in sitelist) {
tdata=mydata[which(mydata$Unique_site_year==s),]
sitesp=unique(tdata$Scientific_name)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$Scientific_name==spp)])
Output[grep(s, Output$Row.names),grep(spp,colnames(Output))] <- abund
}
}

Output[is.na(Output)] <- 0

#### Add flood frequency zones and flood history to data frame


