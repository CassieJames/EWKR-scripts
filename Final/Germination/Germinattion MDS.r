###################################################################################################################################
#### Germination trials

library(vegan)
library(labdsv)
library(ggplot2)
library(gridExtra)
library(indicspecies)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mydata=read.csv("Germination trials.csv") 
spinfo=read.csv("Species_info.csv")
mergeddata=merge(mydata,spinfo,by.x="Species...rectified", by.y="Germination.species.list", all.x=TRUE)
write.csv(mergeddata, file = "Germination merged data checkV2.csv") # save data out 


# Create matrix
species=unique(mydata$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

tada [is.na(tada )] <- 0

# Dealing with species that occurred in the controls

#controls=tada[grep("Control", rownames(tada)), ]
#controls=controls[, colSums(controls)>0]

tada=tada[  -which(rownames(tada) %in% c("LM_the real CC","MM_Control", "Control")), ]
tada=tada[ , -which(colnames(tada) %in% c("Spirodela spp.", "Adiantum sp.", "Cardamine flexuosa")) ]

tada2<-dropspc(tada, 2)   # remove species that occur in = < two sites/date combos
#tada2 <-decostand(tada2,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative values
tadaPA<-decostand(tada2,method="pa")

result<-metaMDS((tadaPA), distance="bray", autotransform=F, k=3,trymax=100)

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3])
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


png(paste("Wetlands together NMDS PA.png",sep=''),width=25, height=20, units='cm', res=500, pointsize=10, bg='white')
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



###################################################################################################################################################################
#### Narran

wetlist=c("NL", "MQ", "MM", "LM")


tada3=tada[grep("NL",rownames(tada)),]
tada4<-dropspc(tada3, 2)   # remove species that occur in = < two sites/date combos
#tada4 <-decostand(tada4,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative values
tadaPA<-decostand(tada4,method="pa")

result<-metaMDS((tadaPA), distance="bray", autotransform=T, k=3,trymax=100)
tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3])

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

NL.FF<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=Flood_Frequency), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$Flood_Frequency, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Narran Lakes, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)		
	
NL.Veg<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= FF)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$Overstory, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Narran Lakes, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)		


###################################################################################################################################################################
#### Marshes

wetlist=c("NL", "MQ", "MM", "LM")


tada3=tada[grep("MQ",rownames(tada)),]

tada4<-dropspc(tada3, 2)   # remove species that occur in = < two sites/date combos
#tada4 <-decostand(tada4,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative values
tadaPA<-decostand(tada4,method="pa")
result<-metaMDS((tadaPA), distance="bray", autotransform=T, k=3,trymax=100)
tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3])

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
	
MQ.FF<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=Flood_Frequency), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$group, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Maquarie Marshes, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)
	
MQ.Veg<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= FF)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$group, size=2)+	
	annotate("text", x = xpos, y=ypos, label = paste("Maquarie Marshes, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)	

###################################################################################################################################################################
#### Mid Murray

wetlist=c("NL", "MQ", "MM", "LM")


tada3=tada[grep("MM",rownames(tada)),]

tada4<-dropspc(tada3, 2)   # remove species that occur in = < two sites/date combos
#tada4 <-decostand(tada4,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative values
tadaPA<-decostand(tada4,method="pa")
result<-metaMDS((tadaPA), distance="bray", autotransform=T, k=3,trymax=100)
tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3])

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
	
MM.FF<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=Flood_Frequency), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$group, size=2)+
    annotate("text", x = xpos,y=ypos, label = paste("Mid Murray, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)			

MM.Veg<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= FF)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$group, size=2)+
    annotate("text", x = xpos, y=ypos, label = paste("Mid Murray, 3D, stress=",round(result$stress,digits = 2)),size=4,hjust=0)		

##############################################################################################################################################################################
#### Lower Murray
wetlist=c("NL", "MQ", "MM", "LM")


tada3=tada[grep("LM",rownames(tada)),]

tada4<-dropspc(tada3, 2)   # remove species that occur in = < two sites/date combos
#tada4 <-decostand(tada4,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative values
tadaPA<-decostand(tada4,method="pa")
result<-metaMDS((tadaPA), distance="bray", autotransform=T, k=2,trymax=100)
tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])

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
	annotate("text", x = xpos, y=ypos, label = paste("Lower Murray, 2D, stress=",round(result$stress,digits = 2)),size=4,hjust = 0)	

LM.FF<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= Overstory)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_FF, aes(x=MDS1, y=MDS2,colour=Flood_Frequency), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.FF$MDS1,y=NMDS.mean.FF$MDS2,label=NMDS.mean.FF$group, size=2)+
	annotate("text",  x = xpos, y=ypos, label = paste("Lower Murray, 2D, stress=",round(result$stress,digits = 2)),size=4,hjust = 0)	

LM.Veg<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape= FF)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Veg, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Veg$MDS1,y=NMDS.mean.Veg$MDS2,label=NMDS.mean.Veg$group, size=2)+
	annotate("text", x = xpos, y=ypos, label = paste("Lower Murray, 2D, stress=",round(result$stress,digits = 2)),size=4,hjust = 0)

	

###############################################################################################################################################################################


png(paste("Wetlands separate NMDS FF PA.png",sep=''),width=25, height=20, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

grid_arrange_shared_legend(NL.FF+nt,MQ.FF+nt,MM.FF+nt,LM.FF+nt, ncol=2,nrow=2)

dev.off()

png(paste("Wetlands separate NMDS Overstory PA.png",sep=''),width=25, height=20, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

grid_arrange_shared_legend(NL.Veg,MQ.Veg,MM.Veg,LM.Veg, ncol=2,nrow=2)

dev.off()

###############################################################################################################################################################################
# Indicator species of wetlands

Wetland=rownames(tada2)
Wetland[grep("NL",Wetland)] <-"Narran Lakes"
Wetland[grep("MQ",Wetland)] <-"Maquarie Marshes"
Wetland[grep("MM",Wetland)] <-"Mid Murray"
Wetland[grep("LM",Wetland)] <-"Lower Murray"
Wetland=as.factor(Wetland)
tada2=as.data.frame(tada2)

indval = multipatt(tada2, Wetland, control = how(nperm=999))
WetlandIndicators <-(summary(indval, alpha=0.05))


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
mydata=read.csv("Germination trials.csv") 
spinfo=read.csv("Species_info.csv")
mergeddata=merge(mydata,spinfo,by.x="Species...rectified", by.y="Germination.species.list", all.x=TRUE)

remove= c("Spirodela spp.", "Adiantum sp.", "Cardamine flexuosa")
mergeddata=mergeddata[-which(mergeddata$Species...rectified %in% remove),]
remove=c("LM_the real CC","MM_Control", "Control")
mergeddata=mergeddata[-which(mergeddata$Site %in% remove),]

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

tada [is.na(tada )] <- 0
tada.all=tada

################################
# Create matrix for exotics
species=unique(mergeddata.exotic$Species...rectified)
sites <- unique(mydata$Site)
tada = matrix(NA,nrow=length(sites),ncol=length(species))#define the output matrix
rownames(tada)=sites
colnames(tada)=species

# Fill matrix
for (s in sites) {

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

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

tdata=mydata[grep(s,mydata$Site),]
sitesp=unique(tdata$Species...rectified)

for (sp in sitesp) {
spdata=tdata[grep(sp,tdata$Species...rectified),]
tada[grep(s,rownames(tada)),grep(sp,colnames(tada))] <-sum(spdata$Count)
}}

tada [is.na(tada )] <- 0
tada.Sedge.rush=tada

###############################

mytada=as.data.frame(rownames(tada))
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


write.csv(mytada, file = "Germination univariate summary.csv") # save data out 



########################################
#### Plots
library(dplyr)
library(ggplot2)
library(gridExtra)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mytada=read.csv("Germination univariate summary.csv") 

variables=c("Abund", "Rich", "ExoticAbund", "ExoticRich","NativeAbund", "NativeRich")

myd=mytada


png(paste('Germination_Seedling Abundances.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(Abund))
 
    ggplot(myd, aes(x = factor(Flow_cat), y = Abund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean seedling abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	

png(paste('Germination_Seedling Richness.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(Rich))

    ggplot(myd, aes(x = factor(Flow_cat), y = Rich, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean species Richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Exotic Abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(ExoticAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = ExoticAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean exotic abundance") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Exotic proportion.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(ExoticAbund/Abund))

    ggplot(myd, aes(x = factor(Flow_cat), y = ExoticAbund/Abund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean exotic proportions") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Exotic Rich.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(ExoticRich))

    ggplot(myd, aes(x = factor(Flow_cat), y = ExoticRich, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean exotic species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Exotic proportion.png',sep=''), width=2000, height=2000, units="px", res=300)	
	mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(Prop.exotic))

    ggplot(myd, aes(x = factor(Flow_cat), y = Prop.exotic, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean proportion of exotics") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Native abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(NativeAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = NativeAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean native abundance") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Native rich.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(NativeRich))

    ggplot(myd, aes(x = factor(Flow_cat), y = NativeRich, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean native species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Annual abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(AnnualAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = AnnualAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean annual abundance") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Annual rich.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(AnnualRich))

    ggplot(myd, aes(x = factor(Flow_cat), y = AnnualRich, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean annual species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Perennial abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(PerennialAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = PerennialAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean perennial abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Forb rich.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(ForbRich))

    ggplot(myd, aes(x = factor(Flow_cat), y = AnnualRich, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean forb species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Forb abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(ForbAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = ForbAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean forb abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
		
png(paste('Germination_Grass rich.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(GrassRich))

    ggplot(myd, aes(x = factor(Flow_cat), y = GrassRich, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean grass species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Grass abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(GrassAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = GrassAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean grass abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_sedge rich.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(Sedge.rush))

    ggplot(myd, aes(x = factor(Flow_cat), y = Sedge.rush, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean sedge species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Sedge abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(Sedge.rushAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = Sedge.rushAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean sedge abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Tree rich.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(TreeRich))

    ggplot(myd, aes(x = factor(Flow_cat), y = TreeRich, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean tree species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Tree abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(TreeAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = TreeAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean tree abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_SS rich.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(SSRich))

    ggplot(myd, aes(x = factor(Flow_cat), y = SSRich, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean subshrub species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_SS abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(SSAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = SSAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean subshurb abundances") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()
	
png(paste('Germination_Shrub rich.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(ShrubRich))

    ggplot(myd, aes(x = factor(Flow_cat), y = ShrubRich, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
    theme_bw()+
    xlab("Flow category")+ # for the x axis label
    ylab("Mean shrub species richness") + labs(colour = "Overstorey")+facet_wrap(. ~ Location,ncol=2)

	dev.off()

png(paste('Germination_Shrub abund.png',sep=''), width=2000, height=2000, units="px", res=300)

mydsum <- ddply(myd,.(Flow_cat,Veg,Location),summarise, val = mean(ShrubAbund))

    ggplot(myd, aes(x = factor(Flow_cat), y = ShrubAbund, colour = Veg)) + 
    geom_point(data = mydsum, aes(y = val)) +
    geom_line(data = mydsum, aes(y = val, group = Veg)) + 
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





