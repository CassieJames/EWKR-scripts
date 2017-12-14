# Scripts to determine antecedent precipitation 
# Written by  C.S.James with code from Tony Ladson (https://tonyladson.wordpress.com/2016/08/15/converting-from-utm-to-lat-long/)
# GNU General Public License .. feel free to use / distribute ... no warranties
# 1st Nov 2016

library(zoo)
library(SDMTools) 
library(dplyr)
source("https://gist.githubusercontent.com/TonyLadson/f37aab3e2ef517188a7f27166307c985/raw/0822970769bc90fcc28052a91b375399d665286e/UTM2deg.R")



data.dir = "/rdsi/vol07/cciaf/Current projects/MD Vegetation/Environmental data"; setwd (data.dir) # set working directory
out.dir="/rdsi/vol07/cciaf/Current projects/MD Vegetation/Environmental data/TERN_data/Precipitation/"
Dates=data.frame(read.csv("HTH_WL_dates_corrected.csv"))
Dates2 <- within(Dates, Date.of.collection <- as.Date((Date.of.collection), format = "%d/%m/%Y")) # make sure r recognises dates as dates
DatesLL <- cbind(Dates2, UTM2deg(Dates2$Easting, Dates$Northing, zone = 55)) # convert UTM to longs and lats


DOI = na.omit(Dates2$Date.of.collection) #create list of dates
DOIUnique=as.Date(unique(na.omit(DOI))) # make sure list is unique
wd = '/rdsi/vol07/ctbcc_data/Climate/CIAS/Australia/1km/baseline.76to05/'
pos = read.csv(paste(wd,'base.positions.csv',sep=""),as.is=TRUE)
data.dir="/rdsi/vol07/ctbcc_data/Climate/eMAST/DAILY/prec/" # location of EMAST tern climate grids
YEARS=c("2008","2009", "2010", "2011", "2012", "2013", "2014") # data only goes to 2012 but will poss be updated soon to 2014
YEARS=c("2013", "2014") # data only goes to 2012 but will poss be updated soon to 2014

for (y in YEARS) {
	
	doi=paste(data.dir, y,sep='')
	foi=list.files(doi)
	
	for(f in foi) { cat(f,'\n') 										# cycle through each file
			tpos=pos
			tasc = read.asc.gz(paste(data.dir,"/",y,"/",f,sep='')) 		# bring in relevant precip asci
			tpos$precip=tasc[cbind(pos$row,pos$col)]     					# Append precip data to temp pos file

			if (f==paste("prec_",y,"0101.asc.gz", sep='')) {
				Jan1=tpos	
			} else {
				Jan1=cbind(Jan1, tpos[,5])
			}
	
	}
	
	jan1back=Jan1	# make temporary copy
    Jan1=round(Jan1,3)	
	tt=gsub("prec_","",foi)
	tt=gsub(".asc.gz","",tt)
	colnames(Jan1) = c('Row',"Col","Lat","Lon", tt)
	write.csv(Jan1,paste(out.dir,y,"_PrecipV2.csv",sep=''),row.names=F) # save as CSV file
	save(Jan1,file=paste(out.dir,y,"_PrecipV2.Rdata",sep='')); rm(Jan1) #write out the data as Rdata object and remove from memory to save space			
	
	}
	

## Read in files and calculate rain history


YEARS=c("2013","2014") # file labels are different for 2013 and 2014 :(

for (y in YEARS) {
	
	doi=paste(data.dir, y,sep='')
	foi=list.files(doi)
	
	for(f in foi) { cat(f,'\n') 										# cycle through each file
			tpos=pos
			tasc = read.asc.gz(paste(data.dir,"/",y,"/",f,sep='')) 		# bring in relevant precip asci
			tpos$precip=tasc[cbind(pos$row,pos$col)]     					# Append precip data to temp pos file

			if (f==paste("ANUClimate_v1-0_rainfall_daily_0-01deg_1970-2014_00000000_",y,"0101.nc.asc.gz", sep='')) {
				Jan1=tpos	
			} else {
				Jan1=cbind(Jan1, tpos[,5])
			}
	
	}
	
	jan1back=Jan1	# make temporary copy
    Jan1=round(Jan1,3)	
	tt=gsub("prec_","",foi)
	tt=gsub(".asc.gz","",tt)
	colnames(Jan1) = c('Row',"Col","Lat","Lon", tt)
	write.csv(Jan1,paste(out.dir,y,"_Precip.csv",sep=''),row.names=F) # save as CSV file
	save(Jan1,file=paste(out.dir,y,"_Precip.Rdata",sep='')); rm(Jan1) #write out the data as Rdata object and remove from memory to save space			
	
	}
	