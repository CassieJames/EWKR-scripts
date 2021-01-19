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
# Multivariate analysis

mydata=tada4 # species abundance metrix 

Wetland=rownames(mydata)
Wetland[grep("NL",Wetland)] <-"Narran Lakes"
Wetland[grep("MQ",Wetland)] <-"Maquarie Marshes"
Wetland[grep("MM",Wetland)] <-"Mid Murray"
Wetland[grep("LM",Wetland)] <-"Lower Murray"
Wetland=as.factor(Wetland)

dis <- vegdist(mydata)
mod <- betadisper(dis, Wetland) # dispersion is signfiicantly lower in MQ but design is unbalanced

# Narran - balanced subset with inland woodlands removed

tada_NL=tada4[grep("NL",rownames(tada4)),]
tada_NL_sub1=tada_NL[-grep("_IW_",rownames(tada_NL)),]

FF_NL=rownames(tada_NL_sub1)

FF_NL[grepl("_C1",FF_NL)] <- "C1"
FF_NL[grepl("_C2",FF_NL)] <- "C2"
FF_NL[grepl("_C3",FF_NL)] <- "C3"
FF_NL[grepl("_C4",FF_NL)] <- "C4"

Veg_NL=rownames(tada_NL_sub1)
Veg_NL[grep("_IS",Veg_NL)] <-"IS"
Veg_NL[grep("_NWW",Veg_NL)] <-"NWW"

Narran.both=cbind(tada_NL_sub1,as.data.frame(FF_NL),as.data.frame(Veg_NL))

Narran.perm <-adonis(tada_NL_sub1 ~ FF_NL*Veg_NL, data=Narran.both ,method = "bray", permutations=999)

dis <- vegdist(tada_NL_sub1)
mod <- betadisper(dis, FF_NL) 

mod <- betadisper(dis, Veg_NL) # F=7.7246, P=0.0084)

# Narran - balanced subset with C2 removed

tada_NL_sub2=tada_NL[-grep("_C2_",rownames(tada_NL)),]

FF_NL=rownames(tada_NL_sub2)

FF_NL[grepl("_C1",FF_NL)] <- "C1"
FF_NL[grepl("_C3",FF_NL)] <- "C3"
FF_NL[grepl("_C4",FF_NL)] <- "C4"

Veg_NL=rownames(tada_NL_sub2)
Veg_NL[grep("_IS",Veg_NL)] <-"IS"
Veg_NL[grep("_NWW",Veg_NL)] <-"NWW"
Veg_NL[grep("_IW",Veg_NL)] <-"IW"

Narran.both=cbind(tada_NL_sub2,as.data.frame(FF_NL),as.data.frame(Veg_NL))

Narran.perm <-adonis(tada_NL_sub2 ~ FF_NL*Veg_NL, data=Narran.both ,method = "bray", permutations=999)

dis <- vegdist(tada_NL_sub2)
mod <- betadisper(dis, FF_NL) 
mod <- betadisper(dis, Veg_NL) 

# Marshes - balanced subset with inland shrubland removed

tada_MQ=tada4[grep("MQ",rownames(tada4)),]
tada_MQ_sub1=tada_MQ[-grep("_IS_",rownames(tada_MQ)),]

#mq.try<-metaMDS(tada_MQ_sub1, distance="bray", autotransform=T, k=2,trymax=100)
#FF_MQ <-as.factor(FF_MQ)
#with(tada_MQ_sub1, ordihull(mq.try, group=FF_MQ))

FF_MQ=rownames(tada_MQ_sub1)

FF_MQ[grepl("_C1",FF_MQ)] <- "C1"
FF_MQ[grepl("_C2",FF_MQ)] <- "C2"
FF_MQ[grepl("_C3",FF_MQ)] <- "C3"
FF_MQ[grepl("_C4",FF_MQ)] <- "C4"

Veg_MQ=rownames(tada_MQ_sub1)
Veg_MQ[grep("_IW",Veg_MQ)] <-"IS"
Veg_MQ[grep("_NWW",Veg_MQ)] <-"NWW"

Marshes.both=cbind(tada_MQ_sub1,as.data.frame(FF_MQ),as.data.frame(Veg_MQ))
Marshes.perm <-adonis(tada_MQ_sub1 ~ FF_MQ*Veg_MQ, data=Marshes.both ,method = "bray", permutations=999)

dis <- vegdist(tada_MQ_sub1)
mod <- betadisper(dis, FF_MQ)

dis <- vegdist(tada_MQ_sub1)
mod <- betadisper(dis, Veg_MQ) 

# Marshes - balanced subset with C1 removed

tada_MQ_sub2=tada_MQ[-grep("_C1_",rownames(tada_MQ)),]

#mq.try<-metaMDS(tada_MQ_sub2, distance="bray", autotransform=T, k=2,trymax=100)
#FF_MQ <-as.factor(FF_MQ)
#with(tada_MQ_sub2, ordihull(mq.try, group=FF_MQ))

FF_MQ=rownames(tada_MQ_sub2)

FF_MQ[grepl("_C2",FF_MQ)] <- "C2"
FF_MQ[grepl("_C3",FF_MQ)] <- "C3"
FF_MQ[grepl("_C4",FF_MQ)] <- "C4"

Veg_MQ=rownames(tada_MQ_sub2)
Veg_MQ[grep("_IW",Veg_MQ)] <-"IW"
Veg_MQ[grep("_NWW",Veg_MQ)] <-"NWW"
Veg_MQ[grep("_IS",Veg_MQ)] <-"IS"

Marshes.both=cbind(tada_MQ_sub2,as.data.frame(FF_MQ),as.data.frame(Veg_MQ))
Marshes.perm <-adonis(tada_MQ_sub2 ~ FF_MQ*Veg_MQ, data=Marshes.both ,method = "bray", permutations=999)

dis <- vegdist(tada_MQ_sub2)
mod <- betadisper(dis, FF_MQ)

TukeyHSD(mod, which = "group", ordered = FALSE, conf.level = 0.95)

dis <- vegdist(tada_MQ_sub2)
mod <- betadisper(dis, Veg_MQ) 

# Middle Murray - balanced subset with C3,C4 and IS removed

tada_MM=tada4[grep("MM",rownames(tada4)),]
tada_MM_sub1=tada_MM[-grep("_C3_",rownames(tada_MM)),]

FF_MM=rownames(tada_MM_sub1)
FF_MM[grepl("_C2",FF_MM)] <- "C2"
FF_MM[grepl("_C1",FF_MM)] <- "C1"

Veg_MM=rownames(tada_MM_sub1)
Veg_MM[grep("_IW",Veg_MM)] <-"IW"
Veg_MM[grep("_NWW",Veg_MM)] <-"NWW"

MM.both=cbind(tada_MM_sub1,as.data.frame(FF_MM),as.data.frame(Veg_MM))
MM.perm <-adonis(tada_MM_sub1 ~ FF_MM*Veg_MM, data=MM.both ,method = "bray", permutations=999)

dis <- vegdist(tada_MM_sub1)
mod <- betadisper(dis, FF_MM)

dis <- vegdist(tada_MM_sub1)
mod <- betadisper(dis, Veg_MM) 

# Lower Murray - balanced subset with C1 removed as absent

tada_LM_sub1=tada4[grep("LM",rownames(tada4)),]


FF_LM=rownames(tada_LM_sub1)
FF_LM[grepl("_C2",FF_LM)] <- "C2"
FF_LM[grepl("_C3",FF_LM)] <- "C3"
FF_LM[grepl("_C4",FF_LM)] <- "C4"

Veg_LM=rownames(tada_LM_sub1)
Veg_LM[grep("_IW",Veg_LM)] <-"IW"
Veg_LM[grep("_NWW",Veg_LM)] <-"NWW"
Veg_LM[grep("_IS",Veg_LM)] <-"IS"

LM.both=cbind(tada_LM_sub1,as.data.frame(FF_LM),as.data.frame(Veg_LM))
LM.perm <-adonis(tada_LM_sub1 ~ FF_LM*Veg_LM, data=LM.both ,method = "bray", permutations=999)

dis <- vegdist(tada_LM_sub1)
mod <- betadisper(dis, FF_LM)

dis <- vegdist(tada_LM_sub1)
mod <- betadisper(dis, Veg_LM) 



 