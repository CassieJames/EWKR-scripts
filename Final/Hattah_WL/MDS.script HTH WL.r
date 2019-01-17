# Script to undertake basic MDS
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016
library(vegan)
library(labdsv)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

#data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.dir="C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
image.dir="C:/Users/jc246980/Documents/MD Vegetation/Plots/"
data.matrix=read.csv("Spp_site_matrix summarised to wetland_HTH_WL_June 2018.csv",row.names = 1 ) # load data - object name is 'tada' ... :)

tada<-dropspc(data.matrix, 2)   # remove species that occur in = < two sites/date combos
tada <-decostand(tada,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative values

result<-metaMDS((tada), distance="bray", autotransform=F, k=3,trymax=100)

#

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$Wetland=rownames(tdata)

tdata$Wetland[grepl("LM_",tdata$Wetland)] <- "LM"
tdata$Wetland[grepl("NL_",tdata$Wetland)] <- "NL"
tdata$Wetland[grepl("MM_",tdata$Wetland)] <- "MM"
tdata$Wetland[grepl("MQ",tdata$Wetland)] <- "MQ"

tdata$Veg<- rownames(tdata)

tdata$Veg[grepl("_IS_",tdata$Veg)] <- "IS"
tdata$Veg[grepl("_NWW_",tdata$Veg)] <- "NWW"
tdata$Veg[grepl("_IW_",tdata$Veg)] <- "IW"

tdata$FF<- rownames(tdata)

tdata$FF[grepl("_C1",tdata$FF)] <- "C1"
tdata$FF[grepl("_C2",tdata$FF)] <- "C2"
tdata$FF[grepl("_C3",tdata$FF)] <- "C3"
tdata$FF[grepl("_C4",tdata$FF)] <- "C4"


tdata$sites[grepl("KT",tdata$sites)] <- "KT"
tdata$sites[grepl("YT",tdata$sites)] <- "YT"
tdata$sites[grepl("MOT",tdata$sites)] <- "MOT"
tdata$sites[grepl("NCT",tdata$sites)] <- "NCT"
tdata$sites[grepl("NN",tdata$sites)] <- "NN"

tdata$group=as.factor(tdata$group)
tdata$sites=as.factor(tdata$sites)


NMDS.mean=aggregate(tdata[,1:2],list(group=tdata$group),mean)
NMDS.mean.sites=aggregate(tdata[,1:2],list(group=tdata$sites),mean)	
	
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





##############################################################
# Code to plot NMDS for each site separately - this code is the 

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT", "KT", "LHAT", "MOT", "NCT", "NN", "YT")



png(paste(image.dir,"Hattah Wetlands NMDS sites separate.png",sep=''),width=20, height=30, units='cm', res=300, pointsize=20, bg='white')
		
ggplot(data = tdata, aes(MDS1, MDS2)) +
geom_point(data=tdata[,c("MDS1","MDS2")],size=1,aes(x=MDS1,y=MDS2),colour="grey")+
geom_point(aes(color = group),size=2)+
facet_wrap(~sites,ncol=2)+geom_path(linetype = "dashed",colour="darkgrey") 


dev.off()


##############################################################
####bioenv


data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Environmental data/Hattah Lakes hydrology/"; setwd (data.dir) 
HydroHTH=data.frame(read.csv("Hydrology_HTH_WL.csv",row.names=1)) # load in hydraulic history for each site

