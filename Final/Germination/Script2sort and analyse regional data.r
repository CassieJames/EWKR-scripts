###################################################################################################################################
#### Germination trials EWKR
#### C James (TropWATER, JCU)
#### No promises!
#### Updated analysis August 2019

library(vegan)
library(labdsv)
library(ggplot2)
library(gridExtra)
library(indicspecies)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mydata=read.csv("Germination trials Aug 2019V2 JN.csv") # this is the dataset that Jason has identified the ambiguous species
mydata=subset(mydata,is.na(mydata$Ambiguous.taxonomy)) # This line subsets the data to only species that are not ambiguous
remove= c("Spirodela spp.", "Adiantum sp.", "Cardamine flexuosa", "No plants")
mydata=subset(mydata,!mydata$Species...rectified %in% remove) 

################################################################################################
# Create first matrix which includes treatments and reps - This is just the basic site by species matrix with everything

species=unique(mydata$Species...rectified)
sites <- unique(mydata$Label) # these are the site names WITH the replicates
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {
tdata=mydata[grep(s,mydata$Label),]
sitesp=unique(tdata$Species...rectified)
for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

tada [is.na(tada )] <- 0

write.csv(tada, file = "Germination species by site.csv") # save data out 

################################################################################################
# Create site by species matrix based on only reps from D1 and S1

dataD1=tada[grep("D1",rownames(tada)),] # takes full matrix (tada) and just subsets to rep 1
dataS1=tada[grep("S1",rownames(tada)),]
dataD1S1=rbind(dataD1,dataS1) # remove data from other replicates as different numbers of replicates undertaken

# Create empty matrix
species=colnames(dataD1S1)
sites <- unique(mydata$Site) # these are the site names WITHOUT the replicates
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix - this adds the results for S1 and d1 together
for (s in sites) {
tdata=dataD1S1[grep(s,rownames(dataD1S1)),] # grab relevant rows from data with only D1 and S1
if(is.null(nrow(tdata))){
tada[grep(s,rownames(tada)),] <-tdata
} else {
tdata2=colSums(tdata) # sum together results from S1 and D1
tada[grep(s,rownames(tada)),] <-tdata2
}}

################################################################################################
# Clean up matrix to remove columns and rows with no records 

tada [is.na(tada )] <- 0
tada2<-dropspc(tada, 2)   # remove species that occur 2 or fewer times
tada3=tada2[rowSums(tada2!= 0) > 0,]
tada4=tada3[,colSums(tada3!= 0) > 0]
tadaPA<-decostand(tada4,method="pa")#  tadaPA on the data with rare species removed otherwise the nMDS is very skewed by a couple of samples

################################################################################################
#### Marshes multivariate

tadaPA_MQ=tadaPA[grep("MQ",rownames(tadaPA)),]  # subset PA to Marshes data only
tadaPA_MQ=tadaPA_MQ[,colSums(tadaPA_MQ!= 0) > 0] # remove species with no data for MQ

tada_MQ=tada4[grep("MQ",rownames(tada4)),]
tada_MQ=tada_MQ[rowSums(tada_MQ!= 0) > 0,]
tada_MQ=tada_MQ[,colSums(tada_MQ!= 0) > 0]

tada_MQ<-dropspc(tada_MQ, 2)   # remove species that occur in = < two sites/date combos within marshes


result<-metaMDS((tada_MQ), distance="bray", autotransform=T, k=3,trymax=100) 
# does not converge very easily so needs a number of runs
tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3]) # extract coordinates for nMDS

#result<-metaMDS((tadaPA_MQ), distance="bray", autotransform=F, k=3,trymax=100)
#tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3])

tdata$Overstory<- rownames(tdata)

tdata$Overstory[grepl("_IS_",tdata$Overstory)] <- "IS"
tdata$Overstory[grepl("_NWW_",tdata$Overstory)] <- "NWW"
tdata$Overstory[grepl("_IW_",tdata$Overstory)] <- "IW"

tdata$FF<- rownames(tdata)

tdata$FF[grepl("_C1",tdata$FF)] <- "C1"
tdata$FF[grepl("_C2",tdata$FF)] <- "C2"
tdata$FF[grepl("_C3",tdata$FF)] <- "C3"
tdata$FF[grepl("_C4",tdata$FF)] <- "C4"

tdata$FFbyVeg <- paste(tdata$FF, "_",tdata$Overstory, sep="")


tdata$Overstory=as.factor(tdata$Overstory)
tdata$FF=as.factor(tdata$FF)
tdata$FFbyVeg=as.factor(tdata$FFbyVeg)

NMDS.mean.FF=aggregate(tdata[,1:2],list(group=tdata$FF),mean)	
NMDS.mean.Veg=aggregate(tdata[,1:2],list(group=tdata$Overstory),mean)	
NMDS.mean.FFbyVeg=aggregate(tdata[,1:2],list(group=tdata$FFbyVeg),mean)	
	
df_ell_Veg <- data.frame()
for(g in levels(tdata$Overstory)){
  df_ell_Veg <- rbind(df_ell_Veg, cbind(as.data.frame(with(tdata[tdata$Overstory==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_FF <- data.frame()
for(g in levels(tdata$FF)){
  df_ell_FF <- rbind(df_ell_FF, cbind(as.data.frame(with(tdata[tdata$FF==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

colnames(df_ell_FF)[3]="Flood_Frequency"

df_ell_FFbyVeg <- data.frame()
for(g in levels(tdata$FFbyVeg)){
  df_ell_FFbyVeg <- rbind(df_ell_FFbyVeg, cbind(as.data.frame(with(tdata[tdata$FFbyVeg==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

ypos = min(tdata$MDS2,na.rm = TRUE) + 0.99*diff(range(tdata$MDS2,na.rm = TRUE))
xpos=min(tdata$MDS1,na.rm = TRUE) + 0.01*diff(range(tdata$MDS1,na.rm = TRUE))
		
MQ<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape = Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FFbyVeg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FFbyVeg$MDS1,y=NMDS.mean.FFbyVeg$MDS2,label=NMDS.mean.FFbyVeg$group, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Maquarie Marshes, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)	
	
MQ.FF<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(colour=FF)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=Flood_Frequency), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$group, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Maquarie Marshes, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)
	
MQ.Veg<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes( colour=Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$group, size=2)	
	#annotate("text", x = xpos, y=ypos, label = paste("3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)	
	
###############################################################################################################################################################################
# Univariate analysis
# For this analysis I do all the calculations to create the summary metrics on full dataset and then subset to site of interest
# This is a bit of a long winded approach as well as I have created separate species by site matrix for each 

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mydata=read.csv("Germination trials Aug 2019V2 JN.csv") # this is the dataset that Jason has identified the ambiguous species
spinfo=read.csv("Species_info_Aug2019.csv")
mergeddata=merge(mydata,spinfo,by.x="Species...rectified", by.y="Germination.species.list", all.x=TRUE) # this merges the germination data with the species info data

remove= c("Spirodela spp.", "Adiantum sp.", "Cardamine flexuosa", "No plants")
mergeddata=subset(mergeddata,!mergeddata$Species...rectified %in% remove) 

dataD1=mergeddata[grep("_D1",mergeddata$Label),]
dataS1=mergeddata[grep("_S1",mergeddata$Label),]
dataD1S1=rbind(dataD1,dataS1) # remove data from other replicates as different numbers of replicates undertaken

mergeddata=dataD1S1 # replaced orginal data with data for only reps 1 - might want to rename rather than override as I have done
mergeddata.wet=subset(mergeddata,mergeddata$Wet.Dry =="Wet")
mergeddata.dry=subset(mergeddata,mergeddata$Wet.Dry =="Dry")

#mergeddata.exotic=mergeddata[(mergeddata$Exotic.Native.x =="Exotic"),]
#mergeddata.native=mergeddata[(mergeddata$Exotic.Native.x =="Native"),]

mergeddata.annual=mergeddata[(mergeddata$Life.history.x =="annual"),]
mergeddata.perennial=mergeddata[(mergeddata$Life.history.x =="perennial"),]

################################
# Create matrix for full species
species=unique(mergeddata$Species...rectified)
sites <- unique(mergeddata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mergeddata[grep(s,mergeddata$Site),c("Species...rectified", "Count")]
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Count")
sitesp=unique(ttdata$Species)

for (sp in sitesp) {

tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-ttdata[ttdata$Species==sp,2]
}
}

tada [is.na(tada )] <- 0
tada.all=tada
tada.all.pa=tada.all
tada.all.pa[tada.all.pa>0] <-1


################################
# Create matrix for annual

species=unique(mergeddata.annual$Species...rectified)
sites <- unique(mergeddata.annual$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
# Fill matrix
for (s in sites) {

tdata=mergeddata.annual[grep(s,mergeddata.annual$Site),c("Species...rectified", "Count")]
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Count")
sitesp=unique(ttdata$Species)

for (sp in sitesp) {

tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-ttdata[ttdata$Species==sp,2]
}
}

tada [is.na(tada )] <- 0
tada.annual=tada
tada.annual.pa=tada.annual
tada.annual.pa[tada.annual.pa>0] <-1

################################
# Create matrix for perennial

species=unique(mergeddata.perennial$Species...rectified)
sites <- unique(mergeddata.perennial$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for (s in sites) {

tdata=mergeddata.perennial[grep(s,mergeddata.perennial$Site),c("Species...rectified", "Count")]
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Count")
sitesp=unique(ttdata$Species)

for (sp in sitesp) {

tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-ttdata[ttdata$Species==sp,2]
}
}
tada [is.na(tada )] <- 0
tada.perennial=tada
tada.perennial.pa=tada.perennial
tada.perennial.pa[tada.perennial.pa>0] <-1

################################
# Create matrix for wet

species=unique(mergeddata.wet$Species...rectified)
sites <- unique(mergeddata.wet$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for (s in sites) {

tdata=mergeddata.wet[grep(s,mergeddata.wet$Site),c("Species...rectified", "Count")]
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Count")
sitesp=unique(ttdata$Species)

for (sp in sitesp) {

tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-ttdata[ttdata$Species==sp,2]
}
}
tada [is.na(tada )] <- 0
tada.wet=tada
tada.wet.pa=tada.wet
tada.wet.pa[tada.wet.pa>0] <-1

################################
# Create matrix for dry

species=unique(mergeddata.dry$Species...rectified)
sites <- unique(mergeddata.dry$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

for (s in sites) {

tdata=mergeddata.dry[grep(s,mergeddata.dry$Site),c("Species...rectified", "Count")]
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Count")
sitesp=unique(ttdata$Species)

for (sp in sitesp) {

tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-ttdata[ttdata$Species==sp,2]
}
}
tada [is.na(tada )] <- 0
tada.dry=tada
tada.dry.pa=tada.dry
tada.dry.pa[tada.dry.pa>0] <-1

###############################
# create new dataframe with summary metrics

mytada=as.data.frame(rownames(tada.all))
mytada$Abund=rowSums(tada.all)
mytada$Rich=rowSums(tada.all.pa)
colnames(mytada)[1]<-"Site"

AnnualAbund=(as.data.frame(rowSums(tada.annual)))
colnames(AnnualAbund)<-c("Abund_annual")
mytada=merge(mytada,AnnualAbund,by.x="Site", by.y="row.names",all.x=TRUE)

AnnualRich=(as.data.frame(rowSums(tada.annual.pa)))
colnames(AnnualRich)<-c("Rich_annual")
mytada=merge(mytada,AnnualRich,by.x="Site", by.y="row.names",all.x=TRUE)

PerennialAbund=(as.data.frame(rowSums(tada.perennial)))
colnames(PerennialAbund)<-c("Abund_perennial")
mytada=merge(mytada,PerennialAbund,by.x="Site", by.y="row.names",all.x=TRUE)

PerennialRich=(as.data.frame(rowSums(tada.perennial.pa)))
colnames(PerennialRich)<-c("Rich_perennial")
mytada=merge(mytada,PerennialRich,by.x="Site", by.y="row.names",all.x=TRUE)

WetAbund=(as.data.frame(rowSums(tada.wet)))
colnames(WetAbund)<-c("Abund_wet")
mytada=merge(mytada,WetAbund,by.x="Site", by.y="row.names",all.x=TRUE)

WetRich=(as.data.frame(rowSums(tada.wet.pa)))
colnames(WetRich)<-c("Rich_wet")
mytada=merge(mytada,WetRich,by.x="Site", by.y="row.names",all.x=TRUE)

DryAbund=(as.data.frame(rowSums(tada.dry)))
colnames(DryAbund)<-c("Abund_dry")
mytada=merge(mytada,DryAbund,by.x="Site", by.y="row.names",all.x=TRUE)

DryRich=(as.data.frame(rowSums(tada.dry.pa)))
colnames(DryRich)<-c("Rich_dry")
mytada=merge(mytada,DryRich,by.x="Site", by.y="row.names",all.x=TRUE)

#### Summary
colnames(mytada)[1]="Site"
mytada$Flow_cat=mytada$Site

mytada$Flow_cat=as.character(mytada$Flow_cat)
mytada$Flow_cat[grep("C1",mytada$Flow_cat)] <- "C1"
mytada$Flow_cat[grepl("C2_",mytada$Flow_cat)] <- "C2"
mytada$Flow_cat[grepl("C3_",mytada$Flow_cat)] <- "C3"
mytada$Flow_cat[grepl("C4_",mytada$Flow_cat)] <- "C4"
mytada$Flow_cat=as.factor(mytada$Flow_cat)

mytada$Location=mytada$Site
mytada$Location=as.character(mytada$Location)
mytada$Location[grep("NL",mytada$Location)] <- "NL"
mytada$Location[grep("MQ",mytada$Location)] <- "MQ"
mytada$Location[grep("MM",mytada$Location)] <- "MM"
mytada$Location[grep("LM",mytada$Location)] <- "LM"
mytada$Location=as.factor(mytada$Location)

mytada$Veg=mytada$Site
mytada$Veg=as.character(mytada$Veg)
mytada$Veg[grep("IS",mytada$Veg)] <- "IS"
mytada$Veg[grep("IW",mytada$Veg)] <- "IW"
mytada$Veg[grep("NWW",mytada$Veg)] <- "NWW"
mytada$Veg=as.factor(mytada$Veg)

mytada[is.na(mytada)] <- 0
write.csv(mytada, file = "Germination univariate summary August 2019.csv") # save data out 

#####################################################################################################
#### Subset summary metrics to site of interest

mytada_MQ=mytada[grep("MQ",mytada$Site),]
myd=mytada_MQ


mydsum.abund <- ddply(myd,.(Flow_cat,Veg),summarise,N=length(Abund), val = mean(Abund), se=sd(Abund)/sqrt(N))
mydsum.abund$metric<-"Abundance"
 
    MQ.abund <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean seedling abundances") +facet_wrap(. ~ Veg,ncol=3)



mydsum.rich <- ddply(myd,.(Flow_cat,Veg),summarise, N=length(Rich), val = mean(Rich), se=sd(Rich)/sqrt(N))
mydsum.rich$metric<-"Richness"

    MQ.rich <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean species Richness") +facet_wrap(. ~ Veg,ncol=3)


mydsum.abund.annual <- ddply(myd,.(Flow_cat,Veg),summarise, N=length(Abund_annual), val = mean(Abund_annual), se=sd(Abund_annual)/sqrt(N))
mydsum.abund.annual$metric<-"Annual abundance"

    MQ.abund.annual <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean annual abundance") + facet_wrap(. ~ Veg,ncol=3)


mydsum.rich.annual <- ddply(myd,.(Flow_cat,Veg),summarise, N=length(Rich_annual), val = mean(Rich_annual), se=sd(Rich_annual)/sqrt(N))
mydsum.rich.annual$metric<-"Annual richness"

    MQ.rich.annual <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean annual species richness") + facet_wrap(. ~ Veg,ncol=3)


mydsum.abund.perennial <- ddply(myd,.(Flow_cat,Veg),summarise, N=length(Abund_perennial), val = mean(Abund_perennial), se=sd(Abund_perennial)/sqrt(N))
mydsum.abund.perennial$metric<-"Perennial abundance"

    MQ.abund.perennial <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean perennial abundances") +facet_wrap(. ~ Veg,ncol=3)


mydsum.rich.perennial <- ddply(myd,.(Flow_cat,Veg),summarise, N=length(Rich_perennial), val = mean(Rich_perennial), se=sd(Rich_perennial)/sqrt(N))
mydsum.rich.perennial$metric<-"Perennial richness"

    MQ.rich.perennial <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean perennial species richness") +facet_wrap(. ~ Veg,ncol=3)


mydsum.abund.wet <- ddply(myd,.(Flow_cat,Veg),summarise, N=length(Abund_wet), val = mean(Abund_wet), se=sd(Abund_wet)/sqrt(N))
mydsum.abund.wet$metric<-"Wet abundance"

    MQ.abund.wet <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean wet species abundances") +facet_wrap(. ~ Veg,ncol=3)
	

mydsum.rich.wet <- ddply(myd,.(Flow_cat,Veg),summarise, N=length(Rich_wet), val = mean(Rich_wet), se=sd(Rich_wet)/sqrt(N))
mydsum.rich.wet$metric<-"Wet richness"

    MQ.rich.wet <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean wet species richness") + facet_wrap(. ~ Veg,ncol=3)


mydsum.abund.dry <- ddply(myd,.(Flow_cat,Veg),summarise, N=length(Abund_dry), val = mean(Abund_dry), se=sd(Abund_dry)/sqrt(N))
mydsum.abund.dry$metric<-"Dry abundance"

    MQ.abund.dry <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean dry abundances") +facet_wrap(. ~ Veg,ncol=3)


mydsum.rich.dry <- ddply(myd,.(Flow_cat,Veg),summarise, N=length(Rich_dry), val = mean(Rich_dry), se=sd(Rich_dry)/sqrt(N))
mydsum.rich.dry$metric<-"Dry richness"

    MQ.rich.dry <-ggplot(mydsum, aes(x = factor(Flow_cat), y = val)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +		
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean dry species richness") +facet_wrap(. ~ Veg,ncol=3)


Facet_plot =rbind(mydsum.abund,mydsum.rich,mydsum.abund.annual,mydsum.rich.annual,mydsum.abund.perennial,mydsum.rich.perennial,mydsum.abund.wet,mydsum.rich.wet,mydsum.abund.dry,mydsum.rich.dry)
Facet_plot$metric=as.factor(Facet_plot$metric)	
Facet_plot$Veg=as.factor(Facet_plot$Veg)	
	
   MQ.univar.metrics <- 
   ggplot(Facet_plot, aes(factor(Flow_cat))) + 
    geom_point(data = Facet_plot, aes(y = val)) +
    geom_line(data = Facet_plot, aes(y = val,group=1)) + 
	geom_errorbar(data=Facet_plot, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ facet_grid(metric ~ Veg,scales="free",switch="y")
	
shortlist= c("Abundance", "Richness", "Annual abundance", "Annual richness","Perennial abundance", "Perennial richness")
Facet_plot1=subset(Facet_plot,Facet_plot$metric %in% shortlist) 
neworder <- c("Richness","Abundance","Annual abundance", "Annual richness","Perennial abundance", "Perennial richness") # force order of plot in facet
Facet_plot1 <- arrange(transform(Facet_plot1,
             metric=factor(metric,levels=neworder)),metric)

MQ.univar.metrics1 <- 
   ggplot(Facet_plot1, aes(factor(Flow_cat))) + 
    geom_point(data = Facet_plot1, aes(y = val)) +
    geom_line(data = Facet_plot1, aes(y = val,group=1)) + 
	geom_errorbar(data=Facet_plot1, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ facet_grid(metric ~ Veg,scales="free",switch="y")
	
shortlist2= c("Wet abundance", "Wet richness","Dry abundance", "Dry richness")
Facet_plot2=subset(Facet_plot,Facet_plot$metric %in% shortlist2) 

MQ.univar.metrics2 <- 
   ggplot(Facet_plot2, aes(factor(Flow_cat))) + 
    geom_point(data = Facet_plot2, aes(y = val)) +
    geom_line(data = Facet_plot2, aes(y = val,group=1)) + 
	geom_errorbar(data=Facet_plot2, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ facet_grid(metric ~ Veg,scales="free",switch="y")

###############################################################################################################################################################################
# different univariate plot arrangement examples
library(gridExtra)

png(paste("Marshes univariate plots.png",sep=''),width=10, height=30, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

grid.arrange(MQ.abund,  MQ.rich, MQ.abund.annual,MQ.rich.annual,MQ.abund.perennial,MQ.rich.perennial, ncol=1,nrow=6)

dev.off()

png(paste("Marshes univariate plots together.png",sep=''),width=20, height=30, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

MQ.univar.metrics


dev.off()
	
	
png(paste("Marshes univariate plots part 1.png",sep=''),width=30, height=30, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

MQ.univar.metrics1


dev.off()
	
		
png(paste("Marshes univariate plots part 2.png",sep=''),width=30, height=30, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

MQ.univar.metrics2


dev.off()

#### THE END