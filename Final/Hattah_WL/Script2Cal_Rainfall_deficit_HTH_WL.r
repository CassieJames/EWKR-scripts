# Scripts to determine antecedent rainfall conditions for sites
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016

library(zoo)
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data"; setwd (data.dir) # set working directory
Dates=data.frame(read.csv("HTH_WL_dates_corrected.csv"))
rain.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/Hattah Lakes rainfall/"; setwd (rain.dir)
Rainfall=data.frame(read.csv("Rainfall_HTH_PARKS.csv"))

Dates2 <- within(Dates, Date.of.collection <- as.Date((Date.of.collection), format = "%d/%m/%Y")) # make sure r recognises dates as dates
Rainfall2 <- within(Rainfall, Date <- as.Date(as.character(Date), format = "%d/%m/%Y"))
Rainfall2[is.na(Rainfall2)] <- 0

DOI = na.omit(Dates2$Date.of.collection) #create list of dates
DOIUnique=as.Date(unique(na.omit(DOI))) # make sure list is unique

Output= matrix(NA,nrow=length(DOIUnique), ncol=1)
colnames(Output)=c("Date")
Output=as.data.frame(Output)
Output$Date<-na.omit(DOIUnique)

# loops don't cope with dates - they convert to numeric so its necessary to use an apply function instead

myfunction<-function (d){

d30 =d-30 	# one month ago
d90=d-90	# three months ago
d180=d-180	# six months ago
d365=d-365	# twelve months ago

loi30=Rainfall2[Rainfall2$Date %in% as.Date(d:d30),]
loi90=Rainfall2[Rainfall2$Date %in% as.Date(d:d90),]
loi180=Rainfall2[Rainfall2$Date %in% as.Date(d:d180),]
loi365=Rainfall2[Rainfall2$Date %in% as.Date(d:d365),]

Rain30=sum(na.omit(loi30$Rain))
Rain90=sum(na.omit(loi90$Rain))
Rain180=sum(na.omit(loi180$Rain))
Rain365=sum(na.omit(loi365$Rain))

Raindef30=sum(na.omit(loi30$No_deficit))
Raindef90=sum(na.omit(loi90$No_deficit))
Raindef180=sum(na.omit(loi180$No_deficit))
Raindef365=sum(na.omit(loi365$No_deficit))
output=c(Rain30, Rain90, Rain180, Rain365,Raindef30,Raindef90,Raindef180,Raindef365 )
return(output)
}

Rainfall.dat=t(sapply(as.list(DOIUnique), myfunction))
Rainfall.dat=as.data.frame(Rainfall.dat)
Rainfall.dat$Date=DOIUnique
colnames(Rainfall.dat)=c("d30", "d90", "d180", "d365","def30", "def90", "def180", "def365", "Date")

write.csv(Rainfall.dat , file = "Rainfall_HTH_WL.csv")