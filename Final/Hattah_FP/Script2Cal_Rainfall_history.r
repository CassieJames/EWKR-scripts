# Scripts to determine antecedent rainfall conditions for sites
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016

library(zoo)
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data"; setwd (data.dir) # set working directory
Dates=data.frame(read.csv("HTH_FP_dates.csv"))
Rainfall=data.frame(read.csv("IDCJAC0009_076043_1800_Data.csv"))

Dates2 <- within(Dates, Sample.date <- as.Date((Sample.date), format = "%d/%m/%Y")) # make sure r recognises dates as dates
Rainfall2 <- within(Rainfall, Date <- as.Date(as.character(Date), format = "%d/%m/%Y"))

DOI = na.omit(Dates2$Sample.date) #create list of dates
DOIUnique=as.Date(unique(na.omit(DOI))) # make sure list is unique

Output= matrix(NA,nrow=length(DOI), ncol=1)
colnames(Output)=c("Date")
Output=as.data.frame(Output)
Output$Date<-na.omit(Dates2$Sample.date)

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

Rain30=sum(na.omit(loi30$Rainfall))
Rain90=sum(na.omit(loi90$Rainfall))
Rain180=sum(na.omit(loi180$Rainfall))
Rain365=sum(na.omit(loi365$Rainfall))
output=c(Rain30, Rain90, Rain180, Rain365)
return(output)
}

paaa=t(sapply(as.list(DOIUnique), myfunction))
paaa=as.data.frame(paaa)
paaa$Date=DOIUnique
colnames(paaa)=c("d30", "d90", "d180", "d365", "Date")

Rainfall.dat <- merge(Output, paaa, by="Date", all.x=TRUE)

write.csv(Rainfall.dat , file = "Rainfall_HTH_FP.csv")