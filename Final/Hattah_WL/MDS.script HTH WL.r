# Script to undertake basic MDS
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016
library(vegan)
library(labdsv)
library(RColorBrewer)

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Spp_site_matrix summarised to wetland_HTH_WL_Dec 2017.csv") # load data - object name is 'tada' ... :)
head(tada)
tada<-dropspc(tada, 2)   # remove species that occur in = < two sites/date combos


result<-metaMDS((tada), distance="bray", autotransform=F, k=2)

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$group=rownames(tdata)

tdata$group[grepl("_08",tdata$group)] <- "Y2008"
tdata$group[grepl("_09",tdata$group)] <- "Y2009"
tdata$group[grepl("_10",tdata$group)] <- "Y2010"
tdata$group[grepl("_11",tdata$group)] <- "Y2011"
tdata$group[grepl("_12",tdata$group)] <- "Y2012"
tdata$group[grepl("_13",tdata$group)] <- "Y2013"
tdata$group[grepl("_14",tdata$group)] <- "Y2014"
tdata$group[grepl("_16",tdata$group)] <- "Y2016"
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

tdata$group=as.factor(tdata$group)
tdata$sites=as.factor(tdata$sites)


NMDS.mean=aggregate(tdata[,1:2],list(group=tdata$group),mean)
				   
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

ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(color = group)) +
    geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)



ord<-ordiellipse(result, tdata$sites, display = "sites", 
                   kind = "se", conf = 0.95, label = T)

NMDS.mean=aggregate(tdata[,1:2],list(group=tdata$sites),mean)
				   
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
				   
df_ell <- data.frame()
for(g in levels(tdata$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(tdata[tdata$group==g,],
                  veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}

ggplot(data = tdata, aes(NMDS1, NMDS2)) + geom_point(aes(color = sites)) +
    geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)









