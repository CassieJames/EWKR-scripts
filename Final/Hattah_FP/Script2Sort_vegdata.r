# Script to turn data into a site by species matrix and to add environmental attributes
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016

# Import data
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv"))

#create matrix for results (species x sites)
specieslist=unique(mydata$Scientific_name)
sitelist=unique(mydata$Unique_site_year)
Output= matrix(NA,nrow=length(unique(mydata$Unique_site_year)), ncol=length(unique(mydata$Scientific_name)))
rownames(Output)=unique(mydata$Unique_site_year)
colnames(Output)=sort(unique(mydata$Scientific_name))


# This section just creates a species x site matrix with no environmental attributes
# for each site and sample date retrieve species list and then loop through and place in results matrix
for(s in sitelist) {
tdata=mydata[which(mydata$Unique_site_year==s),]
sitesp=unique(tdata$Scientific_name)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$Scientific_name==spp)])
Output[grep(s, rownames(Output)),grep(spp,colnames(Output))] <- abund
}
}

# Replace all NAs with zero for analysis
Output[is.na(Output)] <- 0

# This section creates a matrix with the environmental attributes tagged onto the end
# Create species matrix with environmental data attributes
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("HTH_FP.csv")) # load species data
image.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/"; setwd(image.dir)
Rainfall.dat=data.frame(read.csv("Rainfall_HTH_FP.csv")) # load rainfall data which is already sorted into dates x rainfall metrics
Rainfall.dat <- within(Rainfall.dat, Date <- as.Date(as.character(Date), format = "%d/%m/%Y")) # ensure dates are recognised in rainfall data
mydata <- within(mydata, Date.of.collection <- as.Date(as.character(Date.of.collection), format = "%d/%m/%Y")) # ensure dates are recognised in species data
Rainfall.dat=Rainfall.dat[!duplicated(Rainfall.dat), ] # remove duplicates
mydata_env=merge(mydata, Rainfall.dat, by.x="Date.of.collection", by.y="Date") # merge species data with rainfall data by date to create a unique site-sample date by rainfall metrics table
mydata_rainfall=mydata_env[,c("Unique_site_year", "d30", "d90", "d180", "d365")]
mydata_rainfall=unique(mydata_rainfall)

# Create species x site matrix and tag on environmental data
specieslist=unique(mydata$Scientific_name)
sitelist=unique(mydata$Unique_site_year)
Output= matrix(NA,nrow=length(unique(mydata$Unique_site_year)), ncol=length(unique(mydata$Scientific_name)))
rownames(Output)=unique(mydata$Unique_site_year)
colnames(Output)=sort(unique(mydata$Scientific_name))
Output=as.data.frame(Output)
Output=merge(Output, mydata_rainfall, by.x="row.names", by.y="Unique_site_year", all.x=TRUE) # add rainfall metrics to empty site x species matrix

for(s in sitelist) { # Fill matrix
tdata=mydata[which(mydata$Unique_site_year==s),]
sitesp=unique(tdata$Scientific_name)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$Scientific_name==spp)])
Output[grep(s, Output$Row.names),grep(spp,colnames(Output))] <- abund
}
}

Output[is.na(Output)] <- 0

Spp_Env_matrix_HTH_FP=Output # take a copy

write.csv(Spp_Env_matrix , file = "Spp_Env_matrix_HTH_FP.csv") # save data out


