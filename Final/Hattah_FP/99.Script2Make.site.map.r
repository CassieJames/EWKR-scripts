


library(zoo)
library(dplyr)
library(ggmap)
library(gridExtra)
source("https://gist.githubusercontent.com/TonyLadson/f37aab3e2ef517188a7f27166307c985/raw/0822970769bc90fcc28052a91b375399d665286e/UTM2deg.R")

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/"; setwd(data.dir)
Dates=data.frame(read.csv("HTH_FP_dates_locs.csv"))
Dates2 <- within(Dates, Date.of.collection <- as.Date((Date.of.collection), format = "%d/%m/%Y")) # make sure r recognises dates as dates
DatesLL <- cbind(Dates2, UTM2deg(Dates2$Easting, Dates$Northing, zone = 54)) # convert UTM to longs and lats

#map locations to check 
map_center <- DatesLL %>% summarise(lon = mean(lon), lat = mean(lat))
map<-ggmap(get_map(map_center, zoom = 11, source = 'google', maptype = 'satellite'), extent = 'device') +
  geom_point(data = DatesLL, aes(lon, lat), colour = 'red', size = 0.2)

png(paste(data.dir,'Map of site locations for Hattah FP.png',sep=''), width=1500, height=1500, units="px", res=500)
grid.arrange(map, ncol=1, nrow =1)
dev.off()
  