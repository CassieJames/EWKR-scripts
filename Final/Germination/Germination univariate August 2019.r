###############################################################################################################################################################################
# Univariate data

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mydata=read.csv("Germination trials Aug 2019.csv") 
spinfo=read.csv("Species_info_Aug2019.csv")
mergeddata=merge(mydata,spinfo,by.x="Species...rectified", by.y="Germination.species.list", all.x=TRUE)

remove= c("Spirodela spp.", "Adiantum sp.", "Cardamine flexuosa")
mergeddata=mergeddata[!mergeddata$Species...rectified %in% remove,]

dataD1=mergeddata[grep("_D1",mergeddata$Label),]
dataS1=mergeddata[grep("_S1",mergeddata$Label),]
dataD1S1=rbind(dataD1,dataS1) # remove data from other replicates as different numbers of replicates undertaken

dataD1S1=(dataD1S1[,c(1:19,54)])

mergeddata=dataD1S1

write.csv(mergeddata, file = "Germination merged to check Workshop.csv") # save data out 

mergeddata.wet=subset(mergeddata,mergeddata$Wet.Dry =="Wet")
mergeddata.dry=subset(mergeddata,mergeddata$Wet.Dry =="Dry")

mergeddata.exotic=mergeddata[(mergeddata$Exotic.Native.x =="Exotic"),]
mergeddata.native=mergeddata[(mergeddata$Exotic.Native.x =="Native"),]

mergeddata.annual=mergeddata[(mergeddata$Life.history.x =="annual"),]
mergeddata.perennial=mergeddata[(mergeddata$Life.history.x =="perennial"),]

mergeddata.Tree=mergeddata[(mergeddata$Life.form.x =="Tree"),]
mergeddata.SS=mergeddata[(mergeddata$Life.form.x =="Sub-shrub"),]
mergeddata.Shrub=mergeddata[(mergeddata$Life.form.x =="Shrub"),]
mergeddata.Forb=mergeddata[(mergeddata$Life.form.x =="Forb"),]
mergeddata.Grass=mergeddata[(mergeddata$Life.form.x =="Grass"),]
mergeddata.Sedge.rush=mergeddata[(mergeddata$Life.form.x =="Sedge/rush"),]

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
# create new dataframe

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
# univariate analysis

##### NARRAN C1 - C4 by IS and NWW

Narran.dat=subset(mytada,mytada$Location =="NL")
Narran.dat=subset(Narran.dat,Narran.dat$Veg =="IS" |Narran.dat$Veg =="NWW" )

##### Abundance - looks like slight effect of flow class but no interaction regardless of method
myvariable=Narran.dat$Abund

NL.myvariable.glmeFull<- glm.nb(myvariable ~ Flow_cat * Veg, data = Narran.dat)
NL.myvariable.glmeFullNoInt<- glm.nb(myvariable ~ Flow_cat + Veg, data = Narran.dat)
NL.myvariable.glmeFullVeg<- glm.nb(myvariable ~Veg, data = Narran.dat)
lrtest(NL.myvariable.glmeFull, NL.myvariable.glmeFullNoInt)
lrtest(NL.myvariable.glmeNoInt, NL.myvariable.glmeFullVeg)

aovFull=aov(log(myvariable+1) ~ Flow_cat * Veg, data = Narran.dat)
aovnoInt=aov(log(myvariable+1) ~ Flow_cat + Veg, data = Narran.dat)
aovNull=aov(log(myvariable+1) ~ 1, data = Narran.dat)

##### Richness - NO EFFECTS
myvariable=Narran.dat$Rich

NL.myvariable.glmeNull<- glm(myvariable ~ 1, data = Narran.dat, family  = poisson)
NL.myvariable.glmeVeg <- glm(myvariable ~ Veg, data = Narran.dat, family  = poisson)
NL.myvariable.glmeFlow <- glm(myvariable ~ Flow_cat, data = Narran.dat, family  = poisson)
NL.myvariable.glmeNoInt <- glm(myvariable ~ Flow_cat + Veg , data = Narran.dat, family  = poisson)
NL.myvariable.glmeFull<- glm(myvariable ~ Flow_cat * Veg, data = Narran.dat, family  = poisson)
NL.AICs <- AIC(NL.myvariable.glmeNull, NL.myvariable.glmeVeg, NL.myvariable.glmeFlow, NL.myvariable.glmeNoInt, NL.myvariable.glmeFull)
lrtest(NL.myvariable.glmeFull, NL.myvariable.glmeNull) # NO SIGNIFICANT FACTORS

##### Abundance_Annual 
myvariable=Narran.dat$Abund_annual

NL.myvariable.glmeFull<- glm.nb(myvariable ~ Flow_cat * Veg, data = Narran.dat)
NL.myvariable.glmeFullNoInt<- glm.nb(myvariable ~ Flow_cat + Veg, data = Narran.dat)
NL.myvariable.glmeFullVeg<- glm.nb(myvariable ~Veg, data = Narran.dat)
lrtest(NL.myvariable.glmeFull, NL.myvariable.glmeFullNoInt)
lrtest(NL.myvariable.glmeNoInt, NL.myvariable.glmeFullVeg)

##### Abundance_richness 
myvariable=Narran.dat$Rich_annual

NL.myvariable.glmeNull<- glm(myvariable ~ 1, data = Narran.dat, family  = poisson)
NL.myvariable.glmeVeg <- glm(myvariable ~ Veg, data = Narran.dat, family  = poisson)
NL.myvariable.glmeFlow <- glm(myvariable ~ Flow_cat, data = Narran.dat, family  = poisson)
NL.myvariable.glmeNoInt <- glm(myvariable ~ Flow_cat + Veg , data = Narran.dat, family  = poisson)
NL.myvariable.glmeFull<- glm(myvariable ~ Flow_cat * Veg, data = Narran.dat, family  = poisson)
NL.AICs <- AIC(NL.myvariable.glmeNull, NL.myvariable.glmeVeg, NL.myvariable.glmeFlow, NL.myvariable.glmeNoInt, NL.myvariable.glmeFull)
lrtest(NL.myvariable.glmeFull, NL.myvariable.glmeNull) # 
lrtest(NL.myvariable.glmeFull, NL.myvariable.glmeNoInt) # 



























########################################
#### Plots

library(ggplot2)
library(gridExtra)
library(plyr)
library(dplyr)
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mytada=read.csv("Germination univariate summary August 2019.csv") 

variables=c("Abund", "Rich", "ExoticAbund", "ExoticRich","NativeAbund", "NativeRich")

myd=mytada


png(paste('Germination_Seedling Abundances August 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise,N=length(Abund), val = mean(Abund), se=sd(Abund)/sqrt(N))
 
    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean seedling abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	

png(paste('Germination_Seedling Richness August 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Rich), val = mean(Rich), se=sd(Rich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean species Richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

	
png(paste('Germination_Annual abund August 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Abund_annual), val = mean(Abund_annual), se=sd(Abund_annual)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean annual abundance") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Annual rich Aug 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Rich_annual), val = mean(Rich_annual), se=sd(Rich_annual)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean annual species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_perennial abund August 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Abund_perennial), val = mean(Abund_perennial), se=sd(Abund_perennial)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean perennial abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_perennial rich August 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Rich_perennial), val = mean(Rich_perennial), se=sd(Rich_perennial)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean perennial species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_wet abund August 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Abund_wet), val = mean(Abund_wet), se=sd(Abund_wet)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean wet species abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
		
png(paste('Germination_wet rich August 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Rich_wet), val = mean(Rich_wet), se=sd(Rich_wet)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean wet species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_dry abund August 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Abund_dry), val = mean(Abund_dry), se=sd(Abund_dry)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean dry abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_dry rich August 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Rich_dry), val = mean(Rich_dry), se=sd(Rich_dry)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +		
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean dry species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
