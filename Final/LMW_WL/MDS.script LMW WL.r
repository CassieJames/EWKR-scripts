# Script to undertake basic MDS
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016
library(vegan)
library(labdsv)
library(RColorBrewer)

#data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/LMW_data_csvs/"; setwd(data.dir)
data.dir="C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Spp_site_matrix summarised to wetland_LMW_WL_June 2018.csv",row.names=1) # load data - object name is 'Output' ... :)
Output=data.matrix
tada<-dropspc(Output, 2)   # remove species that occur in = < two sites/date combos
tada=tada[rowSums(Output!= 0) > 0,]	# remove sites with no records 
tada <-decostand(tada,method="total", margin=1) # this standardises each site by its total - means that differences amongst sites are based on relative value

result<-metaMDS((tada), distance="bray", autotransform=F, k=3,trymax=100,noshare=0.1,path="shortest")

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2], MDS3=result$points[,3])
tdata$group=rownames(tdata)

tdata$group[grepl("_08",tdata$group)] <- "2008"
tdata$group[grepl("_09",tdata$group)] <- "2009"
tdata$group[grepl("_10",tdata$group)] <- "2010"
tdata$group[grepl("_11",tdata$group)] <- "2011"
tdata$group[grepl("_12",tdata$group)] <- "2012"
tdata$group[grepl("_13",tdata$group)] <- "2013"
tdata$group[grepl("_14",tdata$group)] <- "2014"
tdata$group[grepl("_16",tdata$group)] <- "2016"
tdata$sites<- rownames(tdata)
tdata$sites[grepl("BB",tdata$sites)] <- "BB"
tdata$sites[grepl("CR",tdata$sites)] <- "CR"
tdata$sites[grepl("WW",tdata$sites)] <- "WW"
tdata$sites[grepl("LP",tdata$sites)] <- "LP"
tdata$sites[grepl("MUH",tdata$sites)] <- "MUH"
tdata$sites[grepl("UL",tdata$sites)] <- "UL"
tdata$sites[grepl("UMWC",tdata$sites)] <- "UMWC"
tdata$sites[grepl("MLH",tdata$sites)] <- "MLH"
tdata$sites[grepl("BI",tdata$sites)] <- "BI"
tdata$sites[grepl("W33",tdata$sites)] <- "W33"
tdata$sites[grepl("SCB",tdata$sites)] <- "SCB"
tdata$sites[grepl("WL",tdata$sites)] <- "WL"

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

png(paste(image.dir,"LMW Wetlands NMDS.png",sep=''),width=25, height=12, units='cm', res=500, pointsize=10, bg='white')
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
# Code to plot NMDS for each site separately


png(paste(image.dir,"LMW NMDS sites separate.png",sep=''),width=20, height=30, units='cm', res=300, pointsize=20, bg='white')
		
ggplot(data = tdata, aes(MDS1, MDS2)) +
geom_point(data=tdata[,c("MDS1","MDS2")],size=1,aes(x=MDS1,y=MDS2),colour="grey")+
geom_point(aes(color = group),size=2)+
facet_wrap(~sites,ncol=2)+geom_path(linetype = "dashed",colour="darkgrey") 

dev.off()


