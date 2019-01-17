
#### Mesocosm study

library(vegan)
library(labdsv)
library(ggplot2)
library(gridExtra)
library(indicspecies)
library(lemon)


data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Mesocosm study/"; setwd(data.dir)
mydata=read.csv("AliveV2.csv") 
mydata.alive=mydata[mydata$mortality=="l",]
mydata.P1=mydata[mydata$sampling.code=="P1",]
mydata.P2=mydata[mydata$sampling.code=="P2",]


mydata.est=mydata[mydata.alive$sampling.code=="est",]
rownames(mydata.est)=paste(mydata.est$species, "_", mydata.est$TT_new,"_",mydata.est$Number,sep="")
mydata.est=mydata.est[mydata.est$sampling.code=="est",c(2:9)]

result<-metaMDS((mydata.est), distance="bray", autotransform=T, k=2,trymax=100)

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$Species=rownames(tdata)

tdata$Species[grepl("RRG_",tdata$Species)] <- "RRG"
tdata$Species[grepl("BB_",tdata$Species)] <-"BB"
tdata$Species[grepl("C_",tdata$Species)] <-"C"


tdata$Treatments=rownames(tdata)

tdata$Treatments[grepl("FFFF",tdata$Treatments)] <- "FFFF"
tdata$Treatments[grepl("FDFD",tdata$Treatments)] <- "FDFD"
tdata$Treatments[grepl("FDDD",tdata$Treatments)] <- "FDDD"
tdata$Treatments[grepl("DDDD",tdata$Treatments)] <- "Control"
tdata$Treatments[grepl("DDFD",tdata$Treatments)] <- "DDFD"

tdata$Species=as.factor(tdata$Species)
tdata$Treatments=as.factor(tdata$Treatments)

NMDS.mean.Species=aggregate(tdata[,1:2],list(group=tdata$Species),mean)
NMDS.mean.Treatments=aggregate(tdata[,1:2],list(group=tdata$Treatments),mean)	

df_ell_Species <- data.frame()
for(g in levels(tdata$Species)){
  df_ell_Species <- rbind(df_ell_Species, cbind(as.data.frame(with(tdata[tdata$Species==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_Treatments <- data.frame()
for(g in levels(tdata$Treatments)){
  df_ell_Treatments <- rbind(df_ell_Treatments, cbind(as.data.frame(with(tdata[tdata$Treatments==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}



P1<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape = Species)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Treatments, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Treatments$MDS1,y=NMDS.mean.Treatments$MDS2,label=NMDS.mean.Treatments$group)

################################################################################################################
	
mydata.P1=mydata.alive[mydata.alive$sampling.code=="P1",]
rownames(mydata.P1)=paste(mydata.P1$species, "_", mydata.P1$TT_new,"_",mydata.P1$Number,sep="")
mydata.P1=mydata.P1[mydata.P1$sampling.code=="P1",c(2:9)]

result<-metaMDS((mydata.P1), distance="bray", autotransform=T, k=2,trymax=100)

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$Species=rownames(tdata)

tdata$Species[grepl("RRG_",tdata$Species)] <- "RRG"
tdata$Species[grepl("BB_",tdata$Species)] <-"BB"
tdata$Species[grepl("C_",tdata$Species)] <-"C"


tdata$Treatments=rownames(tdata)

tdata$Treatments[grepl("FFFF",tdata$Treatments)] <- "FFFF"
tdata$Treatments[grepl("FDFD",tdata$Treatments)] <- "FDFD"
tdata$Treatments[grepl("FDDD",tdata$Treatments)] <- "FDDD"
tdata$Treatments[grepl("DDDD",tdata$Treatments)] <- "Control"
tdata$Treatments[grepl("DDFD",tdata$Treatments)] <- "DDFD"

tdata$Species=as.factor(tdata$Species)
tdata$Treatments=as.factor(tdata$Treatments)

NMDS.mean.Species=aggregate(tdata[,1:2],list(group=tdata$Species),mean)
NMDS.mean.Treatments=aggregate(tdata[,1:2],list(group=tdata$Treatments),mean)	

df_ell_Species <- data.frame()
for(g in levels(tdata$Species)){
  df_ell_Species <- rbind(df_ell_Species, cbind(as.data.frame(with(tdata[tdata$Species==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_Treatments <- data.frame()
for(g in levels(tdata$Treatments)){
  df_ell_Treatments <- rbind(df_ell_Treatments, cbind(as.data.frame(with(tdata[tdata$Treatments==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}


P2<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape = Species)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Treatments, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Treatments$MDS1,y=NMDS.mean.Treatments$MDS2,label=NMDS.mean.Treatments$group)
	
vec.sp<-envfit(result$points, mydata.P1, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)

vec.sp.df=vec.sp.df[, c(1:2)]*0.6
vec.sp.df$species<-rownames(vec.sp.df)
	
P2.1<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape = Treatments)) +theme_bw(base_size = 12)+
    geom_path(data=df_ell_Species, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=1)+
	geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
    arrow = arrow(length = unit(0.3, "cm")),colour="grey") + 
    geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=species),size=3)+
    coord_fixed()+ guides(colour=FALSE)+coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
	annotate("text",x=NMDS.mean.Species$MDS1,y=NMDS.mean.Species$MDS2,label=NMDS.mean.Species$group)+
	theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank())

	################################################################################################
	
mydata.P2=mydata.alive[mydata.alive$sampling.code=="P2",]
rownames(mydata.P2)=paste(mydata.P2$species, "_", mydata.P2$TT_new,"_",mydata.P2$Number,sep="")
mydata.P2=mydata.P2[mydata.P2$sampling.code=="P2",c(2:9)]



result<-metaMDS((mydata.P2), distance="bray", autotransform=T, k=2,trymax=100)

tdata=data.frame(MDS1=result$points[,1],MDS2=result$points[,2])
tdata$Species=rownames(tdata)

tdata$Species[grepl("RRG_",tdata$Species)] <- "RRG"
tdata$Species[grepl("BB_",tdata$Species)] <-"BB"
tdata$Species[grepl("C_",tdata$Species)] <-"C"


tdata$Treatments=rownames(tdata)

tdata$Treatments[grepl("FFFF",tdata$Treatments)] <- "FFFF"
tdata$Treatments[grepl("FDFD",tdata$Treatments)] <- "FDFD"
tdata$Treatments[grepl("FDDD",tdata$Treatments)] <- "FDDD"
tdata$Treatments[grepl("DDDD",tdata$Treatments)] <- "Control"
tdata$Treatments[grepl("DDFD",tdata$Treatments)] <- "DDFD"

tdata$Species=as.factor(tdata$Species)
tdata$Treatments=as.factor(tdata$Treatments)

NMDS.mean.Species=aggregate(tdata[,1:2],list(group=tdata$Species),mean)
NMDS.mean.Treatments=aggregate(tdata[,1:2],list(group=tdata$Treatments),mean)	

df_ell_Species <- data.frame()
for(g in levels(tdata$Species)){
  df_ell_Species <- rbind(df_ell_Species, cbind(as.data.frame(with(tdata[tdata$Species==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

df_ell_Treatments <- data.frame()
for(g in levels(tdata$Treatments)){
  df_ell_Treatments <- rbind(df_ell_Treatments, cbind(as.data.frame(with(tdata[tdata$Treatments==g,],
                  veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}


P3<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape = Species)) +theme_classic(base_size = 12)+
    geom_path(data=df_ell_Treatments, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=2)+
    annotate("text",x=NMDS.mean.Treatments$MDS1,y=NMDS.mean.Treatments$MDS2,label=NMDS.mean.Treatments$group)

	
vec.sp<-envfit(result$points, mydata.P2, perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)

vec.sp.df=vec.sp.df[, c(1:2)]*0.6
vec.sp.df$species<-rownames(vec.sp.df)
	
P3.1<-ggplot(data = tdata, aes(MDS1, MDS2)) + geom_point(aes(shape = Treatments)) +theme_bw(base_size = 12)+
    geom_path(data=df_ell_Species, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=1)+
	geom_segment(data=vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
    arrow = arrow(length = unit(0.3, "cm")),colour="grey") + 
    geom_text(data=vec.sp.df,aes(x=MDS1,y=MDS2,label=species),size=3)+	coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) +
    guides(colour=FALSE)+
	annotate("text",x=NMDS.mean.Species$MDS1,y=NMDS.mean.Species$MDS2,label=NMDS.mean.Species$group)+
	theme(plot.background = element_blank() ,panel.grid.major = element_blank(),panel.grid.minor = element_blank())
	
png(paste("Mesocosm nMDS plots V1.png",sep=''),width=25, height=14, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(2,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns
	
	grid_arrange_shared_legend(P2.1, P3.1, ncol = 2)
	
	dev.off()

	
	
	
	
	
	
	
	
	
	