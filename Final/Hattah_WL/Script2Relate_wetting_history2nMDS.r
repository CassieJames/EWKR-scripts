# Script to relate wetting history to species composition
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# June 2018
###########################################################################################
library(vegan)
library(labdsv)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(dplyr)

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
image.dir="C:/Users/jc246980/Documents/MD Vegetation/Plots/"
Output=read.csv("Spp_site_year_transect matrix HTH_WL_June 2018.csv",row.names = 1 ) 

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")

#data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd(data.dir)
data.dir="C:/Users/jc246980/Documents/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd(data.dir)
Hydrodata=data.frame(read.csv("Hydraulics_HTH_WL.csv",row.names=1)) # load in hydraulic history for each site


Hydrodata$wetyear <-"FALSE"
Hydrodata$wetyear[which(Hydrodata$WetYear>0)]<-"TRUE"


for(w in wetlist) { # 

if(w=="BIT"){
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="BRT"){
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="BOT"){
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="BLT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
}
if(w=="CCS"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2011",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2012",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}

if(w=="NCT"){
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
}
if(w=="HT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2012",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="KT"){
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2014",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
}
if(w=="LHAT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
}
if(w=="MOT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="NN"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2010",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2011",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2012",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2014",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
if(w=="YT"){
Hydrodata[grepl("2008",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2009",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2011",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2012",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2013",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-FALSE
Hydrodata[grepl("2014",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
Hydrodata[grepl("2016",Hydrodata$HydroYear)& grepl(w,Hydrodata$Site.ID) ,c("wetyear")]<-TRUE
}
}



Hydrodata$TSLWcats<-">1000"
Hydrodata$TSLWcats[which(Hydrodata$TSLW<=90)]<-"<=90"
Hydrodata$TSLWcats[which(Hydrodata$TSLW>90)]<-">90"
Hydrodata$TSLWcats[which(Hydrodata$TSLW>182)]<-">182"
Hydrodata$TSLWcats[which(Hydrodata$TSLW>365)]<-">365"
Hydrodata$TSLWcats[which(Hydrodata$TSLW>730)]<-">730"

Hydrodata$TSLWcats=as.factor(Hydrodata$TSLWcats)


######################################################################################################
## Script to determine cv of maximum depth

Species_hydro=merge(Hydrodata,Output, by.x="Unique_site_year", by.y="row.names") # merge 

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
write.csv(Output , file = "Hydrology_HTH_WL.csv")

Species_hydro$sites<- as.character(Species_hydro$Site.ID)


Species_hydro$sites[grepl("BIT",Species_hydro$sites)] <- "BIT"
Species_hydro$sites[grepl("BLT",Species_hydro$sites)] <- "BLT"
Species_hydro$sites[grepl("BOT",Species_hydro$sites)] <- "BOT"
Species_hydro$sites[grepl("BRT",Species_hydro$sites)] <- "BRT"
Species_hydro$sites[grepl("CCS",Species_hydro$sites)] <- "CCS"
Species_hydro$sites[grepl("HT",Species_hydro$sites)] <- "HT"
Species_hydro$sites[grepl("LHAT",Species_hydro$sites)] <- "LHAT"
Species_hydro$sites[grepl("KT",Species_hydro$sites)] <- "KT"
Species_hydro$sites[grepl("YT",Species_hydro$sites)] <- "YT"
Species_hydro$sites[grepl("MOT",Species_hydro$sites)] <- "MOT"
Species_hydro$sites[grepl("NCT",Species_hydro$sites)] <- "NCT"
Species_hydro$sites[grepl("NN",Species_hydro$sites)] <- "NN"
Species_hydro$HydroYear=as.factor(Species_hydro$HydroYear)


Speciesdat =select(Species_hydro,Abu.theo:Zygo.sp.)
mine<-aggregate(Speciesdat,by = list(Species_hydro$sites,Species_hydro$HydroYear,Species_hydro$wetyear),FUN=sum)
rownames(mine)=paste(mine[,1],"_",mine[,2],"_",mine[,3],sep="")

tada<-dropspc(mine[,-c(1:3)], 2)   # remove species that occur in = < two sites/date combos
tada <-decostand(tada,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative values
tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records 


result<-metaMDS((tada), distance="bray", autotransform=F, k=3,trymax=20, noshare=0.1)

#

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2], MDS3=result$points[,3])
tdata$group=rownames(tdata)

tdata$group[grepl("2008",tdata$group)] <- "2008"
tdata$group[grepl("2009",tdata$group)] <- "2009"
tdata$group[grepl("2010",tdata$group)] <- "2010"
tdata$group[grepl("2011",tdata$group)] <- "2011"
tdata$group[grepl("2012",tdata$group)] <- "2012"
tdata$group[grepl("2013",tdata$group)] <- "2013"
tdata$group[grepl("2014",tdata$group)] <- "2014"
tdata$group[grepl("2015",tdata$group)] <- "2015"
tdata$group[grepl("2016",tdata$group)] <- "2016"
tdata$sites<- rownames(tdata)
tdata$sites[grepl("BIT",tdata$sites)] <- "BIT"
tdata$sites[grepl("BLT",tdata$sites)] <- "BLT"
tdata$sites[grepl("BOT",tdata$sites)] <- "BOT"
tdata$sites[grepl("BRT",tdata$sites)] <- "BRT"
tdata$sites[grepl("CCS",tdata$sites)] <- "CCS"
tdata$sites[grepl("HT",tdata$sites)] <- "HT"
tdata$sites[grepl("LHAT",tdata$sites)] <- "LHAT"
tdata$sites[grepl("KT",tdata$sites)] <- "KT"
tdata$sites[grepl("YT",tdata$sites)] <- "YT"
tdata$sites[grepl("MOT",tdata$sites)] <- "MOT"
tdata$sites[grepl("NCT",tdata$sites)] <- "NCT"
tdata$sites[grepl("NN",tdata$sites)] <- "NN"

tdata$wet<- rownames(tdata)
tdata$wet[grepl("TRUE",tdata$wet)] <- "WET"
tdata$wet[grepl("FALSE",tdata$wet)] <- "DRY"

tdata$group=as.factor(tdata$group)
tdata$sites=as.factor(tdata$sites)
tdata$wet=as.factor(tdata$wet)


NMDS.mean=aggregate(tdata[,1:2],list(group=tdata$group),mean)
NMDS.mean.sites=aggregate(tdata[,1:2],list(group=tdata$sites),mean)	
NMDS.mean.wet=aggregate(tdata[,1:2],list(group=tdata$wet),mean)	
	
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
				   
df_ell <- data.frame()
for(g in levels(tdata$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(tdata[tdata$group==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_sites <- data.frame()
for(g in levels(tdata$sites)){
  df_ell_sites <- rbind(df_ell_sites, cbind(as.data.frame(with(tdata[tdata$sites==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_wet <- data.frame()
for(g in levels(tdata$wet)){
  df_ell_wet <- rbind(df_ell_wet, cbind(as.data.frame(with(tdata[tdata$wet==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}




png(paste(image.dir,"Hattah Wetlands NMDS axes 1 and 2.png",sep=''),width=25, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

yearnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)

sitenmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = sites)) +
    geom_path(data=df_ell_sites, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.sites$MDS1,y=NMDS.mean.sites$MDS2,label=NMDS.mean.sites$group)
	
	grid.arrange(yearnmds,sitenmds, ncol = 2)
	
dev.off()


png(paste(image.dir,"Hattah Wetlands NMDS axes 1 and 2 wet versus dry.png",sep=''),width=12.5, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

wetnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = wet)) +
    geom_path(data=df_ell_wet, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.wet$MDS1,y=NMDS.mean.wet$MDS2,label=NMDS.mean.wet$group)
wetnmds
	
dev.off()

#######################################################################################################################
#### Alternative analysis using TSLW categories


Speciesdat =select(Species_hydro,Abu.theo:Zygo.sp.)
mine<-aggregate(Speciesdat,by = list(Species_hydro$sites,Species_hydro$HydroYear,Species_hydro$TSLWcats),FUN=sum)
rownames(mine)=paste(mine[,1],"_",mine[,2],"_",mine[,3],sep="")

tada<-dropspc(mine[,-c(1:3)], 2)   # remove species that occur in = < two sites/date combos
tada <-decostand(tada,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative values
tada=tada[rowSums(tada!= 0) > 0,]	# remove sites with no records 


result<-metaMDS((tada), distance="bray", autotransform=F, k=3,trymax=20,noshare=0.1) # should be able to use above nMDS


tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2],MDS2=result$points[,3])
tdata$group=rownames(tdata2)

tdata$group[grepl("2008",tdata$group)] <- "2008"
tdata$group[grepl("2009",tdata$group)] <- "2009"
tdata$group[grepl("2010",tdata$group)] <- "2010"
tdata$group[grepl("2011",tdata$group)] <- "2011"
tdata$group[grepl("2012",tdata$group)] <- "2012"
tdata$group[grepl("2013",tdata$group)] <- "2013"
tdata$group[grepl("2014",tdata$group)] <- "2014"
tdata$group[grepl("2015",tdata$group)] <- "2015"
tdata$group[grepl("2016",tdata$group)] <- "2016"
tdata$sites<- rownames(tdata)
tdata$sites[grepl("BIT",tdata$sites)] <- "BIT"
tdata$sites[grepl("BLT",tdata$sites)] <- "BLT"
tdata$sites[grepl("BOT",tdata$sites)] <- "BOT"
tdata$sites[grepl("BRT",tdata$sites)] <- "BRT"
tdata$sites[grepl("CCS",tdata$sites)] <- "CCS"
tdata$sites[grepl("HT",tdata$sites)] <- "HT"
tdata$sites[grepl("LHAT",tdata$sites)] <- "LHAT"
tdata$sites[grepl("KT",tdata$sites)] <- "KT"
tdata$sites[grepl("YT",tdata$sites)] <- "YT"
tdata$sites[grepl("MOT",tdata$sites)] <- "MOT"
tdata$sites[grepl("NCT",tdata$sites)] <- "NCT"
tdata$sites[grepl("NN",tdata$sites)] <- "NN"

tdata$TSLWcats<- rownames(tdata)
tdata$TSLWcats[grepl("<=90",tdata$TSLWcats)] <- "<3 months"
tdata$TSLWcats[grepl(">90",tdata$TSLWcats)] <- ">3 months"
tdata$TSLWcats[grepl(">182",tdata$TSLWcats)] <- ">6 months"
tdata$TSLWcats[grepl(">365",tdata$TSLWcats)] <- ">12 months"
tdata$TSLWcats[grepl(">730",tdata$TSLWcats)] <- ">2 Years"
tdata$TSLWcats[grepl(">1000",tdata$TSLWcats)] <- ">3 Years"

tdata$group=as.factor(tdata$group)
tdata$sites=as.factor(tdata$sites)
tdata$TSLWcats=as.factor(tdata$TSLWcats)


NMDS.mean=aggregate(tdata[,1:2],list(group=tdata$group),mean)
NMDS.mean.sites=aggregate(tdata[,1:2],list(group=tdata$sites),mean)	
NMDS.mean.TSLWcats=aggregate(tdata[,1:2],list(group=tdata$TSLWcats),mean)	
	
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
				   
df_ell <- data.frame()
for(g in levels(tdata$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(tdata[tdata$group==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_sites <- data.frame()
for(g in levels(tdata$sites)){
  df_ell_sites <- rbind(df_ell_sites, cbind(as.data.frame(with(tdata[tdata$sites==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_TSLWcats <- data.frame()
for(g in levels(tdata$TSLWcats)){
  df_ell_TSLWcats <- rbind(df_ell_TSLWcats, cbind(as.data.frame(with(tdata[tdata$TSLWcats==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}




png(paste(image.dir,"Hattah Wetlands NMDS axes 1 and 2 TSLWcats.png",sep=''),width=25, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

yearnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)

sitenmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = sites)) +
    geom_path(data=df_ell_sites, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.sites$MDS1,y=NMDS.mean.sites$MDS2,label=NMDS.mean.sites$group)
	
	grid.arrange(yearnmds,sitenmds, ncol = 2)
	
dev.off()


png(paste(image.dir,"Hattah Wetlands NMDS axes 1 and 2 TSLWcats.png",sep=''),width=12.5, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

TSLWcatsnmds <-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = TSLWcats)) +
    geom_path(data=df_ell_TSLWcats, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.TSLWcats$MDS1,y=NMDS.mean.TSLWcats$MDS2,label=NMDS.mean.TSLWcats$group)
TSLWcatsnmds
	
dev.off()


