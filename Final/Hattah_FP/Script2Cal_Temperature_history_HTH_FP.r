# Scripts to determine temperature conditions
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 6 September 2016


library(zoo)
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data"; setwd (data.dir) # set working directory
Dates=data.frame(read.csv("HTH_FP_dates_corrected.csv"))
data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes temperature/"; setwd (data.dir) # set working directory
Temp=data.frame(read.csv("IDCJAC0010_076047_1800_Data.csv"))

Dates2 <- within(Dates, Date.of.collection <- as.Date((Date.of.collection), format = "%d/%m/%Y")) # make sure r recognises dates as dates
Temp2 <- within(Temp, Date <- as.Date(as.character(Date), format = "%d/%m/%Y"))

DOI = na.omit(Dates2$Date.of.collection) #create list of dates
DOIUnique=as.Date(unique(na.omit(DOI))) # make sure list is unique

Output= matrix(NA,nrow=length(DOI), ncol=1)
colnames(Output)=c("Date")
Output=as.data.frame(Output)
Output$Date<-na.omit(Dates2$Date.of.collection)


myfunction<-function (d){

d30 =d-30 	# one month ago
d90=d-90	# three months ago
d180=d-180	# six months ago
d365=d-365	# twelve months ago


loi30=Temp2[Temp2$Date %in% (d:d30),]
loi90=Temp2[Temp2$Date %in% (d:d90),]
loi180=Temp2[Temp2$Date %in% (d:d180),]
loi365=Temp2[Temp2$Date %in% (d:d365),]

MeanTemp30=mean(na.omit(loi30$Maximum.temperature..Degree.C.))
MaxTemp30=max(na.omit(loi30$Maximum.temperature..Degree.C.))
MinTemp30=min(na.omit(loi30$Maximum.temperature..Degree.C.))

MeanTemp90=mean(na.omit(loi90$Maximum.temperature..Degree.C.))
MaxTemp90=max(na.omit(loi90$Maximum.temperature..Degree.C.))
MinTemp90=min(na.omit(loi90$Maximum.temperature..Degree.C.))

MeanTemp180=mean(na.omit(loi180$Maximum.temperature..Degree.C.))
MaxTemp180=max(na.omit(loi180$Maximum.temperature..Degree.C.))
MinTemp180=min(na.omit(loi180$Maximum.temperature..Degree.C.))

MeanTemp365=mean(na.omit(loi365$Maximum.temperature..Degree.C.))
MaxTemp365=max(na.omit(loi365$Maximum.temperature..Degree.C.))
MinTemp365=min(na.omit(loi365$Maximum.temperature..Degree.C.))

output=c(MeanTemp30, MaxTemp30, MinTemp30,MeanTemp90, MaxTemp90, MinTemp90,MeanTemp180, MaxTemp180, MinTemp180,MeanTemp365, MaxTemp365, MinTemp365)
return(output)
}

paaa=t(sapply(as.list(DOIUnique), myfunction))
paaa=as.data.frame(paaa)
paaa$Date=DOIUnique
colnames(paaa)=c('MeanTemp30', 'MaxTemp30', 'MinTemp30','MeanTemp90', 'MaxTemp90', 'MinTemp90','MeanTemp180', 'MaxTemp180', 'MinTemp180','MeanTemp365', 'MaxTemp365', 'MinTemp365', 'Date')

Temperature.dat <- merge(Output, paaa, by="Date", all.x=TRUE)

write.csv(Temperature.dat , file = "Temperature_HTH_FP.csv")