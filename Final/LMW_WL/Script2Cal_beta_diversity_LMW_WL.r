### Script to determine species accumulation curves and beta diversity calculations for each wetland
### Cassie James
###  7th June 2018

# Import data
#data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 
data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 
mydata=data.frame(read.csv("LMW_WL.csv"))
species=unique(mydata$Scientific.name)
mycodes=data.frame(read.csv("LMW_WL_sp_codes.csv")) # this list has now been amended
mydata=merge(mydata,mycodes,by="Scientific.name")

#Spp_no_codes<-setdiff(unique(mydata$Scientific.name), unique(mydata2$Scientific.name)) # checks for differences to make sure all the species have codes
 
specieslist=unique(mydata$sp_code)
sitelist=unique(mydata$Site_year)
takeout=c("Inundate", "Lea.litt", "Bar.grou") # remove non species
specieslist=specieslist[!specieslist %in% takeout]
Output= matrix(NA,nrow=length(unique(mydata$Site_year)), ncol=length(specieslist))
rownames(Output)=sitelist
colnames(Output)=sort(specieslist)
Output=as.data.frame(Output)

for(s in sitelist) { # Fill matrix
tdata=mydata[which(mydata$Site_year==s),]
sitesp=NULL
abund=NULL
sitesp=unique(tdata$sp_code)
takeout=c("Inundate", "Lea.litt", "Bar.grou")
sitesp=sitesp[!sitesp %in% takeout]
if(length(sitesp)>0){
s=as.character(s)
for (spp in sitesp) {
abund=max(tdata$Abundance[which(tdata$sp_code==spp)])
Output[grep(s, rownames(Output), fixed=TRUE),grep(spp,colnames(Output), fixed=TRUE)] <- abund # needed fixed=TRUE to work
}
}
}
Output[is.na(Output)] <- 0 # replace nas with zeros in species matrix


#### beta diversity within each wetland between years - output indicates very high turnover but thats probably because of a mixture of wet and dry years

wetlist=c("BB","CR","LP","UL","MUH","BI","UMWC","W33","SCB","MLH","WL","WW")


tada = matrix(NA,nrow=length(wetlist),ncol=5)#define the output matrix
rownames(tada)=wetlist
colnames(tada)=c("n","BSor", "BSim", "BNes","BTbc")

right = function(text, num_char) {
  substr(text, nchar(text) - (num_char-1), nchar(text))
}
		
for(w in wetlist) { # 
tdata=Output[grep(w,rownames(Output)),]

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

groups<-rownames(tdata)
groups=paste(right(groups,2),sep="")
groups=as.factor(groups)

mydata <- aggregate(tdata, by = list(groups), FUN = sum)
mean(vegdist(mydata[,-1]))
Baselga<-nestedbetasor(mydata[,-1])
tada[grep(w,rownames(tada)),1] = nrow(mydata)
tada[grep(w,rownames(tada)),2] = Baselga[[3]]
tada[grep(w,rownames(tada)),3] = Baselga[[1]]
tada[grep(w,rownames(tada)),4] = Baselga[[2]]
tada[grep(w,rownames(tada)),5] = mean(vegdist(mydata[,-1]))
}

write.csv(tada , file = "LMW_betawithinsites.csv") # save data out

