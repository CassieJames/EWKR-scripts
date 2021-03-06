###################################################################################################################################
#### Germination trials
#### Updated analysis August 2019

library(vegan)
library(labdsv)
library(ggplot2)
library(gridExtra)
library(indicspecies)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mydata=read.csv("Germination trials Aug 2019V2 JN.csv") 
mydata=subset(mydata,is.na(mydata$Ambiguous.taxonomy))
remove= c("Spirodela spp.", "Adiantum sp.", "Cardamine flexuosa", "No plants")
mydata=subset(mydata,!mydata$Species...rectified %in% remove) 

################################################################################################
# Create first matrix which includes treatments and reps
species=unique(mydata$Species...rectified)
sites <- unique(mydata$Label)
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
# Create matrix based on only reps from D1 and S1

dataD1=tada[grep("D1",rownames(tada)),]
dataS1=tada[grep("S1",rownames(tada)),]
dataD1S1=rbind(dataD1,dataS1) # remove data from other replicates as different numbers of replicates undertaken

#write.csv(dataD1S1, file = "Germination species by site rep1 only.csv") # save data out 

species=colnames(dataD1S1)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix - this adds the results for S1 and d1 together
for (s in sites) {
tdata=dataD1S1[grep(s,rownames(dataD1S1)),] # grab relevant rows
if(is.null(nrow(tdata))){
tada[grep(s,rownames(tada)),] <-tdata
} else {
tdata2=colSums(tdata) # sum together results from S1 and D1
tada[grep(s,rownames(tada)),] <-tdata2
}}

################################################################################################
# Clean up matrix to remove columns and rows with no records

tada [is.na(tada )] <- 0
tada2<-dropspc(tada, 2)   # remove species that occur in = < two sites/date combos
tada3=tada2[rowSums(tada2!= 0) > 0,]
tada4=tada3[,colSums(tada3!= 0) > 0]

tada=tada[rowSums(tada!= 0) > 0,]
tada=tada[,colSums(tada!= 0) > 0]
tadaPA<-decostand(tada4,method="pa") #  tadaPA on the data with rare species removed otherwise the nMDS is very skewed by a couple of samples


################################################################################################
# Run nMDS using PA data
#write.csv(tadaPA, file = "Germination species by site matric rep1 PA.csv") # save data out 

result<-metaMDS((tadaPA), distance="bray", autotransform=F, k=2,trymax=100)

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$Wetland=rownames(tdata)

tdata$Wetland[grepl("LM_",tdata$Wetland)] <- "LM"
tdata$Wetland[grepl("NL_",tdata$Wetland)] <- "NL"
tdata$Wetland[grepl("MM_",tdata$Wetland)] <- "MM"
tdata$Wetland[grepl("MQ",tdata$Wetland)] <- "MQ"

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

tdata$Wetland=as.factor(tdata$Wetland)
tdata$Overstory=as.factor(tdata$Overstory)
tdata$FF=as.factor(tdata$FF)
tdata$FFbyVeg=as.factor(tdata$FFbyVeg)


NMDS.mean.Wetland=aggregate(tdata[,1:2],list(group=tdata$Wetland),mean)
NMDS.mean.FF=aggregate(tdata[,1:2],list(group=tdata$FF),mean)	
colnames(NMDS.mean.FF)[1]="Flood_Frequency"
NMDS.mean.Veg=aggregate(tdata[,1:2],list(group=tdata$Overstory),mean)	
colnames(NMDS.mean.Veg)[1]="Overstory"
NMDS.mean.FFbyVeg=aggregate(tdata[,1:2],list(group=tdata$FFbyVeg),mean)	
	
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
				   
df_ell <- data.frame()
for(g in levels(tdata$Wetland)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(tdata[tdata$Wetland==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}


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

df_ell_FFbyVeg <- data.frame()
for(g in levels(tdata$FFbyVeg)){
  df_ell_FFbyVeg <- rbind(df_ell_FFbyVeg, cbind(as.data.frame(with(tdata[tdata$FFbyVeg==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}


png(paste("Wetlands together NMDS PA August 2019 Species sorted.png",sep=''),width=25, height=20, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

wetlandsnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = Wetland)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Wetland$MDS1,y=NMDS.mean.Wetland$MDS2,label=NMDS.mean.Wetland$group)
	
vegnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = Overstory))  +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$Overstory)
	
FFnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = FF))  +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$Flood_Frequency)	
	
FFbyVegnmds <-ggplot(data = tdata, aes(MDS1, MDS2))  +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FFbyVeg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FFbyVeg$MDS1,y=NMDS.mean.FFbyVeg$MDS2,label=NMDS.mean.FFbyVeg$group, size=2)	
	
	grid.arrange(wetlandsnmds,vegnmds,FFnmds,FFbyVegnmds, ncol = 2)
	
dev.off()


ypos = min(tdata$MDS2,na.rm = TRUE) + 0.99*diff(range(tdata$MDS2,na.rm = TRUE))
xpos=min(tdata$MDS1,na.rm = TRUE) + 0.01*diff(range(tdata$MDS1,na.rm = TRUE))

wetlandsnmds_PA <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = Wetland)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Wetland$MDS1,y=NMDS.mean.Wetland$MDS2,label=NMDS.mean.Wetland$group)+
	annotate("text", x = xpos, y=ypos, label = paste("Presence/Absence, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)	


################################################################################################
# Run nMDS using Abundance data

#tada2 <-decostand(tada2,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative values - 
#write.csv(tada4, file = "Germination species by site matric rep1 Abund.csv") # save data out 


result<-metaMDS((tada4), distance="bray", autotransform=T, k=3,trymax=100)

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$Wetland=rownames(tdata)

tdata$Wetland[grepl("LM_",tdata$Wetland)] <- "LM"
tdata$Wetland[grepl("NL_",tdata$Wetland)] <- "NL"
tdata$Wetland[grepl("MM_",tdata$Wetland)] <- "MM"
tdata$Wetland[grepl("MQ",tdata$Wetland)] <- "MQ"

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

tdata$Wetland=as.factor(tdata$Wetland)
tdata$Overstory=as.factor(tdata$Overstory)
tdata$FF=as.factor(tdata$FF)
tdata$FFbyVeg=as.factor(tdata$FFbyVeg)


NMDS.mean.Wetland=aggregate(tdata[,1:2],list(group=tdata$Wetland),mean)
NMDS.mean.FF=aggregate(tdata[,1:2],list(group=tdata$FF),mean)	
colnames(NMDS.mean.FF)[1]="Flood_Frequency"
NMDS.mean.Veg=aggregate(tdata[,1:2],list(group=tdata$Overstory),mean)	
colnames(NMDS.mean.Veg)[1]="Overstory"
NMDS.mean.FFbyVeg=aggregate(tdata[,1:2],list(group=tdata$FFbyVeg),mean)	
	
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
				   
df_ell <- data.frame()
for(g in levels(tdata$Wetland)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(tdata[tdata$Wetland==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}


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

df_ell_FFbyVeg <- data.frame()
for(g in levels(tdata$FFbyVeg)){
  df_ell_FFbyVeg <- rbind(df_ell_FFbyVeg, cbind(as.data.frame(with(tdata[tdata$FFbyVeg==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}


png(paste("Wetlands together NMDS Abundances August 2019 Species sorted.png",sep=''),width=25, height=20, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

wetlandsnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = Wetland)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Wetland$MDS1,y=NMDS.mean.Wetland$MDS2,label=NMDS.mean.Wetland$group)
	
vegnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = Overstory))  +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$Overstory)
	
FFnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = FF))  +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$Flood_Frequency)	
	
FFbyVegnmds <-ggplot(data = tdata, aes(MDS1, MDS2))  +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FFbyVeg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FFbyVeg$MDS1,y=NMDS.mean.FFbyVeg$MDS2,label=NMDS.mean.FFbyVeg$group, size=2)	
	
	grid.arrange(wetlandsnmds,vegnmds,FFnmds,FFbyVegnmds, ncol = 2)
	
dev.off()

ypos = min(tdata$MDS2,na.rm = TRUE) + 0.99*diff(range(tdata$MDS2,na.rm = TRUE))
xpos=min(tdata$MDS1,na.rm = TRUE) + 0.01*diff(range(tdata$MDS1,na.rm = TRUE))

wetlandsnmds_abund <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = Wetland)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Wetland$MDS1,y=NMDS.mean.Wetland$MDS2,label=NMDS.mean.Wetland$group)+
	annotate("text", x = xpos, y=ypos, label = paste("Abundance, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)	
	
	
	
###################################################################################################################################################################
# Plots of PA and abundance together

library(plyr)
library(lemon)

png(paste("Wetlands NMDS Aug 2019 Abundance and PA.png",sep=''),width=25, height=12.5, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

grid_arrange_shared_legend(wetlandsnmds_PA,wetlandsnmds_abund,ncol=2,nrow=1)

dev.off()


###################################################################################################################################################################
#### Narran

tadaPA_NL=tadaPA[grep("NL",rownames(tadaPA)),]
tadaPA_NL=tadaPA_NL[rowSums(tadaPA_NL!= 0) > 0,]
tadaPA_NL=tadaPA_NL[,colSums(tadaPA_NL!= 0) > 0]

tada_NL=tada4[grep("NL",rownames(tada4)),]
tada_NL=tada_NL[rowSums(tada_NL!= 0) > 0,]
tada_NL=tada_NL[,colSums(tada_NL!= 0) > 0]

tada_NL<-dropspc(tada_NL, 2)   # remove species that occur in = < two sites/date combos

result<-metaMDS((tada_NL), distance="bray", autotransform=T, k=3,trymax=100)
#result<-metaMDS((tadaPA_NL), distance="bray", autotransform=F, k=2,trymax=100)

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3])
#tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])

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
colnames(NMDS.mean.FF)[1]="Flood_Frequency"
NMDS.mean.Veg=aggregate(tdata[,1:2],list(group=tdata$Overstory),mean)	
colnames(NMDS.mean.Veg)[1]="Overstory"
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
		
NL<-ggplot(data = tdata, aes(MDS1, MDS2))+ geom_point(aes(shape= Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FFbyVeg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FFbyVeg$MDS1,y=NMDS.mean.FFbyVeg$MDS2,label=NMDS.mean.FFbyVeg$group, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Narran Lakes, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)		

NL.FF<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes( colour=FF)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=Flood_Frequency), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$Flood_Frequency, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Narran Lakes, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)		
	
NL.Veg<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(colour=Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$Overstory, size=2)
	#annotate("text", x = xpos, y=ypos, label = paste("3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)		

NL.FF.Abund=NL.FF
NL.Veg.Abund=NL.Veg
###################################################################################################################################################################
#### Marshes

tadaPA_MQ=tadaPA[grep("MQ",rownames(tadaPA)),]
tadaPA_MQ=tadaPA_MQ[rowSums(tadaPA_MQ!= 0) > 0,]
tadaPA_MQ=tadaPA_MQ[,colSums(tadaPA_MQ!= 0) > 0]

tada_MQ=tada4[grep("MQ",rownames(tada4)),]
tada_MQ=tada_MQ[rowSums(tada_MQ!= 0) > 0,]
tada_MQ=tada_MQ[,colSums(tada_MQ!= 0) > 0]

tada_MQ<-dropspc(tada_MQ, 2)   # remove species that occur in = < two sites/date combos


result<-metaMDS((tada_MQ), distance="bray", autotransform=T, k=3,trymax=100)
tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3])

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

MQ.FF.abund=MQ.FF
MQ.Veg.abund=MQ.FF
	
###################################################################################################################################################################
#### Mid Murray

tadaPA_MM=tadaPA[grep("MM",rownames(tadaPA)),]
tadaPA_MM=tadaPA_MM[rowSums(tadaPA_MM!= 0) > 0,]
tadaPA_MM=tadaPA_MM[,colSums(tadaPA_MM!= 0) > 0]

tada_MM=tada4[grep("MM",rownames(tada4)),]
tada_MM=tada_MM[rowSums(tada_MM!= 0) > 0,]
tada_MM=tada_MM[,colSums(tada_MM!= 0) > 0]

tada_MM<-dropspc(tada_MM, 2)   # remove species that occur in = < two sites/date combos

result<-metaMDS((tada_MM), distance="bray", autotransform=T, k=2,trymax=100)
tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])

#result<-metaMDS((tadaPA_MM), distance="bray", autotransform=F, k=2,trymax=100)
#tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])

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
		
MM<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FFbyVeg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FFbyVeg$MDS1,y=NMDS.mean.FFbyVeg$MDS2,label=NMDS.mean.FFbyVeg$group, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Mid Murray, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)		

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)
cols=cols[2:4]
	
MM.FF<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(colour=FF)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=Flood_Frequency), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$group, size=2)+
    annotate("text", x = xpos,y=ypos, label = paste("Mid Murray, 2D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)			

MM.Veg<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(colour=Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+scale_color_manual(values=cols)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$group, size=2)
    #annotate("text", x = xpos, y=ypos, label = paste("2D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)	

MM.FF.abund=MM.FF
MM.Veg.abund=MM.Veg	

##############################################################################################################################################################################
#### Lower Murray
tadaPA_LM=tadaPA[grep("LM",rownames(tadaPA)),]
tadaPA_LM=tadaPA_LM[rowSums(tadaPA_LM!= 0) > 0,]
tadaPA_LM=tadaPA_LM[,colSums(tadaPA_LM!= 0) > 0]

tada_LM=tada4[grep("LM",rownames(tada4)),]
tada_LM=tada_LM[rowSums(tada_LM!= 0) > 0,]
tada_LM=tada_LM[,colSums(tada_LM!= 0) > 0]

tada_LM<-dropspc(tada_LM, 2)   # remove species that occur in = < two sites/date combos

result<-metaMDS((tada_LM), distance="bray", autotransform=T, k=3,trymax=100)
tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3])

#result<-metaMDS((tadaPA_LM), distance="bray", autotransform=F, k=3,trymax=100)
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
		
LM<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$group, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Lower Murray, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust = 0)	

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)
cols=cols[2:4]
	
LM.FF<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(colour=FF))+theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=Flood_Frequency), size=1, linetype=2)+scale_color_manual(values=cols)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$group, size=2)+
	annotate("text",  x = xpos, y=ypos, label = paste("Lower Murray, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust = 0)	

LM.Veg<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes( colour=Overstory))+theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$group, size=2)
	#annotate("text", x = xpos, y=ypos, label = paste("2D, stress=",round(result$stress,digits = 2)),size=4,hjust = 0)

LM.FF.abund=LM.FF	
LM.Veg.abund=LM.Veg


###############################################################################################################################################################################
library(plyr)
library(lemon)

png(paste("Wetlands separate NMDS FF August 2019.png",sep=''),width=10, height=30, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

grid_arrange_shared_legend(NL.FF,MQ.FF,MM.FF,LM.FF, ncol=1,nrow=4)

dev.off()

png(paste("Wetlands separate NMDS Overstory August 2019.png",sep=''),width=10, height=30, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

grid_arrange_shared_legend(NL.Veg,MQ.Veg,MM.Veg,LM.Veg, ncol=1,nrow=4)

dev.off()



###############################################################################################################################################################################
# Indicator species of wetlands

Wetbystr= substr(rownames(tadaPA),1,nchar(rownames(tadaPA))-2)

agg = aggregate(tadaPA,
                by = list(Wetbystr),
                FUN = sum)
tadaPAagg<-decostand(agg[,-1],method="pa")
rownames(tadaPAagg)<-agg[,1]
				
Wetland=rownames(tadaPAagg)
Wetland[grep("NL",Wetland)] <-"Narran Lakes"
Wetland[grep("MQ",Wetland)] <-"Maquarie Marshes"
Wetland[grep("MM",Wetland)] <-"Mid Murray"
Wetland[grep("LM",Wetland)] <-"Lower Murray"
Wetland=as.factor(Wetland)
tadaPAagg=as.data.frame(tadaPAagg)

indval = multipatt(tadaPAagg, Wetland, control = how(nperm=999))
WetlandIndicators <-summary(indval, alpha=0.05, indvalcomp=TRUE)


# Indicator species of overstory

Overstory=rownames(tada2)
Overstory[grep("_IS",Overstory)] <-"IS"
Overstory[grep("_IW",Overstory)] <-"IW"
Overstory[grep("_NWW",Overstory)] <-"NWW"

Overstory=as.factor(Overstory)
tada2=as.data.frame(tada2)

indval = multipatt(tada2, Overstory, control = how(nperm=999))
WetlandIndicators <-(summary(indval, alpha=0.05))


###############################################################################################################################################################################
# Univariate data

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mydata=read.csv("Germination trials Aug 2019.csv") 
#spinfo=read.csv("Species_info.csv")
#mergeddata=merge(mydata,spinfo,by.x="Species...rectified", by.y="Germination.species.list", all.x=TRUE)

remove= c("Spirodela spp.", "Adiantum sp.", "Cardamine flexuosa")
#mergeddata=mergeddata[-which(mergeddata$Species...rectified %in% remove),]

dataD1=mydata[grep("_D1",mydata$Label),]
dataS1=mydata[grep("_S1",mydata$Label),]
dataD1S1=rbind(dataD1,dataS1) # remove data from other replicates as different numbers of replicates undertaken
dataD1S1=dataD1S1[!dataD1S1$Species...rectified %in% remove, ]

mergeddata=dataD1S1
mergeddata.exotic=mergeddata[(mergeddata$Exotic.Native =="Exotic"),]
mergeddata.native=mergeddata[(mergeddata$Exotic.Native =="Native"),]

mergeddata.annual=mergeddata[(mergeddata$Life.history =="annual"),]
mergeddata.perennial=mergeddata[(mergeddata$Life.history =="perennial"),]

mergeddata.Tree=mergeddata[(mergeddata$Life.form =="Tree"),]
mergeddata.SS=mergeddata[(mergeddata$Life.form =="Sub-shrub"),]
mergeddata.Shrub=mergeddata[(mergeddata$Life.form =="Shrub"),]
mergeddata.Forb=mergeddata[(mergeddata$Life.form =="Forb"),]
mergeddata.Grass=mergeddata[(mergeddata$Life.form =="Grass"),]
mergeddata.Sedge.rush=mergeddata[(mergeddata$Life.form =="Sedge/rush"),]

################################
# Create matrix for full species
species=unique(mergeddata$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.all=tada

################################
# Create matrix for full species
species=unique(mergeddata.exotic$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.exotic=tada
###############################
# Create matrix for natives

species=unique(mergeddata.native$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.native=tada

################################
# Create matrix for annual

species=unique(mergeddata.annual$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.annual=tada

################################
# Create matrix for perennial

species=unique(mergeddata.perennial$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.perennial=tada

################################
# Create matrix for tree

species=unique(mergeddata.Tree$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.tree=tada
################################
# Create matrix for Shrub

species=unique(mergeddata.Shrub$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.shrub=tada

################################
# Create matrix for SS

species=unique(mergeddata.SS$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.SS=tada

################################
# Create matrix for Forb

species=unique(mergeddata.Forb$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0

tada.Forb=tada

################################
# Create matrix for Grass

species=unique(mergeddata.Grass$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.Grass=tada

################################
# Create matrix for Sedge.rush

species=unique(mergeddata.Sedge.rush$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),c("Site_Treat_Rep","Treatment.rep","Species...rectified", "Count")]
repcount=unique(tdata$Treatment.rep) # return number of replicates for Ds and Ss together
Ds=length(grep("D", repcount))
Ss=length(grep("S", repcount))
ttdata=aggregate(tdata$Count, by=list(tdata$Species...rectified,tdata$Site_Treat_Rep),FUN=sum) # merge records for same species from same replicate
colnames(ttdata)=c("Species", "Site_Treat_Rep", "Count")
tdataD=ttdata[grep("_D",ttdata$Site_Treat_Rep),]
tdataS=ttdata[grep("_S",ttdata$Site_Treat_Rep),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {

spdataD=tdataD[grep(sp,tdataD$Species),]
spdataS=tdataS[grep(sp,tdataS$Species),]

Averaged=(sum(spdataD$Count)/Ds)+(sum(spdataS$Count)/Ss)
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-Averaged
}
}

tada [is.na(tada )] <- 0
tada.Sedge.rush=tada

###############################

mytada=as.data.frame(rownames(tada.all))
mytada$Abund=rowSums(tada.all)
tada.all[tada.all > 0] <- 1 
mytada$Rich=rowSums(tada.all)

mytada$ExoticAbund=rowSums(tada.exotic)
tada.exotic[tada.exotic > 0] <- 1 
mytada$ExoticRich=rowSums(tada.exotic)

mytada$NativeAbund=rowSums(tada.native)
tada.native[tada.native > 0] <- 1 
mytada$NativeRich=rowSums(tada.native)

mytada$AnnualAbund=rowSums(tada.annual)
tada.annual[tada.annual > 0] <- 1 
mytada$AnnualRich=rowSums(tada.annual)

mytada$PerennialAbund=rowSums(tada.perennial)
tada.perennial[tada.perennial > 0] <- 1 
mytada$PerennialRich=rowSums(tada.perennial)

mytada$TreeAbund=rowSums(tada.tree)
tada.tree[tada.tree > 0] <- 1 
mytada$TreeRich=rowSums(tada.tree)

mytada$ShrubAbund=rowSums(tada.shrub)
tada.shrub[tada.shrub > 0] <- 1 
mytada$ShrubRich=rowSums(tada.shrub)

mytada$SSAbund=rowSums(tada.SS)
tada.SS[tada.SS > 0] <- 1 
mytada$SSRich=rowSums(tada.SS)

mytada$ForbAbund=rowSums(tada.Forb)
tada.Forb[tada.Forb > 0] <- 1 
mytada$ForbRich=rowSums(tada.Forb)

mytada$GrassAbund=rowSums(tada.Grass)
tada.Grass[tada.Grass > 0] <- 1 
mytada$GrassRich=rowSums(tada.Grass)

mytada$Sedge.rushAbund=rowSums(tada.Sedge.rush)
tada.Sedge.rush[tada.Sedge.rush > 0] <- 1 
mytada$Sedge.rush=rowSums(tada.Sedge.rush)

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


write.csv(mytada, file = "Germination univariate summary May 2019.csv") # save data out 



########################################
#### Plots

library(ggplot2)
library(gridExtra)
library(plyr)
library(dplyr)
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mytada=read.csv("Germination univariate summary May 2019.csv") 

variables=c("Abund", "Rich", "ExoticAbund", "ExoticRich","NativeAbund", "NativeRich")

myd=mytada


png(paste('Germination_Seedling Abundances May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise,N=length(Abund), val = mean(Abund), se=sd(Abund)/sqrt(N))
 
    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean seedling abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	

png(paste('Germination_Seedling Richness May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Rich), val = mean(Rich), se=sd(Rich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean species Richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Exotic Abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise,  N=length(ExoticAbund), val = mean(ExoticAbund), se=sd(ExoticAbund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val,colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean exotic abundance") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Exotic proportion May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(ExoticAbund/Abund), val = mean(ExoticAbund/Abund), se=sd(ExoticAbund/Abund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean exotic proportions") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Exotic Rich May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(ExoticRich), val = mean(ExoticRich), se=sd(ExoticRich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean exotic species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Native abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise,  N=length(NativeAbund), val = mean(NativeAbund), se=sd(NativeAbund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean native abundance") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Native rich May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(NativeRich), val = mean(NativeRich), se=sd(NativeRich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean native species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Annual abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(AnnualAbund), val = mean(AnnualAbund), se=sd(AnnualAbund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean annual abundance") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Annual rich May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(AnnualRich), val = mean(AnnualRich), se=sd(AnnualRich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean annual species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Perennial abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(PerennialAbund), val = mean(PerennialAbund), se=sd(PerennialAbund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean perennial abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Forb rich May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(ForbRich), val = mean(ForbRich), se=sd(ForbRich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean forb species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Forb abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(ForbAbund), val = mean(ForbAbund), se=sd(ForbAbund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean forb abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
		
png(paste('Germination_Grass rich May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(GrassRich), val = mean(GrassRich), se=sd(GrassRich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean grass species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Grass abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(GrassRich), val = mean(GrassRich), se=sd(GrassRich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean grass abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_sedge rich May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Sedge.rush), val = mean(Sedge.rush), se=sd(Sedge.rush)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +		
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean sedge species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Sedge abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(Sedge.rushAbund), val = mean(Sedge.rushAbund), se=sd(Sedge.rushAbund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +		
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean sedge abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Tree rich May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(TreeRich), val = mean(TreeRich), se=sd(TreeRich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) +
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean tree species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Tree abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise,  N=length(TreeAbund), val = mean(TreeAbund), se=sd(TreeAbund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean tree abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_SS rich May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(SSRich), val = mean(SSRich), se=sd(SSRich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean subshrub species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_SS abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(SSAbund), val = mean(SSAbund), se=sd(SSAbund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean subshurb abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Shrub rich May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(ShrubRich), val = mean(ShrubRich), se=sd(ShrubRich)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean shrub species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Shrub abund May 2019.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, N=length(ShrubAbund), val = mean(ShrubAbund), se=sd(ShrubAbund)/sqrt(N))

    ggplot(mydsum, aes(x = factor(Flow_cat), y = val, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
	geom_errorbar(data=mydsum, aes(ymin=val-se, ymax=val+se), width=.1) +	
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean shrub abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

#########################################################################################################################
#### Subset data into the balnaced component. subset 1 is FF C2, C3 and C4 for Veg IS and NWW at locations LM, MQ and NL


mytada=read.csv("Germination univariate summary.csv") 

subset1 <- mytada[!mytada$Location=="Middle Murray",]
subset1 <- subset1[!subset1$Flow_cat=="C1",]
subset1 <- subset1[!subset1$Veg=="IW",]





