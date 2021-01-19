### Script to determine beta diversity calculations for each wetland
### Cassie James
### 8 June 2018

library(betapart)
library(vegan)

#############################################################################################
### Code to determine beta diversity for all wetland complexes together
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.HLWL=read.csv("Spp_site matrix HTH_WL_May 2019.csv",row.names=1) # load data - object name is 'tada' ... :)
Output.HLWL =data.matrix.HLWL
Output.HLWL=Output.HLWL[rowSums(Output.HLWL!= 0) > 0,]	# remove columns and or rows with no records
Output.HLWL=Output.HLWL[,colSums(Output.HLWL!= 0) > 0]

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/LMW_data_csvs/"; setwd(data.dir)
data.matrix.LMW=read.csv("Spp_site matrix LMW_WL_May 2019.csv", row.names=1) # load data - object name is 'Output' ... :)
Output.LMW =data.matrix.LMW
Output.LMW=Output.LMW[rowSums(Output.LMW!= 0) > 0,]	# remove columns and or rows with no records
Output.LMW=Output.LMW[,colSums(Output.LMW!= 0) > 0]

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
data.matrix.Gun=read.csv("Spp_site_matrix Gunbower_WL_May 2019.csv",row.names=1) # load data - object name is 'tada' ... :)
Output.Gun =data.matrix.Gun
Output.Gun=Output.Gun[rowSums(Output.Gun!= 0) > 0,]	# remove columns and or rows with no records
Output.Gun=Output.Gun[,colSums(Output.Gun!= 0) > 0]


data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/KP_data_csvs/"; setwd (data.dir) 
data.matrix.KP=read.csv("Spp_site_matrix KP_WL_May 2019.csv",row.names=1) # load data - object name is 'tada' ... :)
Output.KP =data.matrix.KP
Output.KP=Output.KP[,-grep("NA",colnames(Output.KP))] # remove 'NA' column
Output.KP=Output.KP[rowSums(Output.KP!= 0) > 0,]	# remove columns and or rows with no records
Output.KP=Output.KP[,colSums(Output.KP!= 0) > 0]

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir) 
data.matrix.Barm=read.csv("Spp_site_matrix Barmah_WL_May 2019.csv",row.names=1) # load data - object name is 'tada' ... :)
Output.Barm =data.matrix.Barm
Output.Barm=Output.Barm[rowSums(Output.Barm!= 0) > 0,]	# remove columns and or rows with no records
Output.Barm=Output.Barm[,colSums(Output.Barm!= 0) >0]

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Chowilla_data_csvs/"; setwd(data.dir)
data.matrix.Chow=read.csv("Spp_site_matrix Chowilla_WL_May 2019.csv",row.names=1) # load data - object name is 'tada' ... :)
Output.Chow =data.matrix.Chow
Output.Chow=Output.Chow[rowSums(Output.Chow!= 0) > 0,]	# remove columns and or rows with no records
Output.Chow=Output.Chow[,colSums(Output.Chow!= 0) >0]
Output.Chow=Output.Chow[,-c(5)] # remove 'NA' column

HLWL.PA<-decostand(Output.HLWL,method="pa")
HLWL.PA.core <- betapart.core(HLWL.PA)
HLWL.PA.multi <- beta.multi(HLWL.PA.core)
HLWL.samp <- beta.sample(HLWL.PA.core, sites=8, samples=1000)
HLWL.dist.n <- HLWL.samp$sampled.values

LMW.PA<-decostand(Output.LMW,method="pa")
LMW.PA.core <- betapart.core(LMW.PA)
LMW.PA.multi <- beta.multi(LMW.PA.core)
LMW.samp <- beta.sample(LMW.PA.core, sites=8, samples=1000)
LMW.dist.n <- LMW.samp$sampled.values

Gun.PA<-decostand(Output.Gun,method="pa")
Gun.PA.core <- betapart.core(Gun.PA)
Gun.PA.multi <- beta.multi(Gun.PA.core)
Gun.samp <- beta.sample(Gun.PA.core, sites=8, samples=1000)
Gun.dist.n <- Gun.samp$sampled.values

KP.PA<-decostand(Output.KP,method="pa")
KP.PA.core <- betapart.core(KP.PA)
KP.PA.multi <- beta.multi(KP.PA.core)
KP.samp <- beta.sample(KP.PA.core, sites=8, samples=1000)
KP.dist.n <- KP.samp$sampled.values

Barm.PA<-decostand(Output.Barm,method="pa")
Barm.PA.core <- betapart.core(Barm.PA)
Barm.PA.multi <- beta.multi(Barm.PA.core)
Barm.samp <- beta.sample(Barm.PA.core, sites=8, samples=1000)
Barm.dist.n <- Barm.samp$sampled.values

Chow.PA<-decostand(Output.Chow,method="pa")
Chow.PA.core <- betapart.core(Chow.PA)
Chow.PA.multi <- beta.multi(Chow.PA.core)
Chow.samp <- beta.sample(Chow.PA.core, sites=8, samples=1000)
Chow.dist.n <- Chow.samp$sampled.values

# colours

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 6
cols = gg_color_hue(n)

png(paste(image.dir,"Wetlands together betapart May 2009.png",sep=''),width=24, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(1,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

plot(density(HLWL.dist.n$beta.SOR), xlim=c(0,0.8), ylim=c(0, 100), xlab="Beta diversity", main="", lwd=3) 

lines(density(HLWL.dist.n$beta.SOR), col=cols[4], lwd=3)
lines(density(HLWL.dist.n$beta.SNE), col=cols[4], lty=1, lwd=2)
lines(density(HLWL.dist.n$beta.SIM), col=cols[4], lty=2, lwd=2)

lines(density(LMW.dist.n$beta.SOR),col=cols[6], lwd=3)
lines(density(LMW.dist.n$beta.SNE),col=cols[6],lty=1, lwd=2)
lines(density(LMW.dist.n$beta.SIM),col=cols[6],lty=2, lwd=2)

lines(density(Gun.dist.n$beta.SOR),col=cols[3], lwd=3)
lines(density(Gun.dist.n$beta.SNE),col=cols[3], lty=1, lwd=2)
lines(density(Gun.dist.n$beta.SIM),col=cols[3], lty=2, lwd=2)

lines(density(KP.dist.n$beta.SOR),col=cols[5], lwd=3)
lines(density(KP.dist.n$beta.SNE),col=cols[5], lty=1, lwd=2)
lines(density(KP.dist.n$beta.SIM),col=cols[5], lty=2, lwd=2)

lines(density(Barm.dist.n$beta.SOR),col=cols[1], lwd=3)
lines(density(Barm.dist.n$beta.SNE),col=cols[1], lty=1, lwd=2)
lines(density(Barm.dist.n$beta.SIM),col=cols[1], lty=2, lwd=2)

lines(density(Chow.dist.n$beta.SOR),col=cols[2], lwd=3)
lines(density(Chow.dist.n$beta.SNE),col=cols[2], lty=1, lwd=2)
lines(density(Chow.dist.n$beta.SIM),col=cols[2], lty=2, lwd=2)

legend(0, 100, legend=c("Barmah", "Chowilla","Gunbower","Hattah", "KP", "LMW"),col=cols, lty=1, lwd=2,cex=1)

dev.off()

##############################################################################################################
#### Code to look at Hattah lakes wetlands differences in beta diversity between wet and dry years

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
mydat=read.csv("Hattah wetlands response by metrics.csv") # load hydrology data
myenv=mydat[,c("Unique_site_year_season","TSLW")]

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
myveg=read.csv("Spp_site_year_season matrix HTH_WL_May 2019 for beta analysis.csv") # save data out

mydata2=merge(myenv, myveg, by.x="Unique_site_year_season",by.y="X",all.y=TRUE) # 

TSLW_1_day <-mydata2[which(mydata2$TSLW<=1),]
TSLW_1_Month <-mydata2[which(mydata2$TSLW<=30 & mydata2$TSLW>1),]
TSLW_3_Month <-mydata2[which(mydata2$TSLW<=90 & mydata2$TSLW>30),]
TSLW_6_Month <-mydata2[which(mydata2$TSLW<=90*2 & mydata2$TSLW>90),]
TSLW_9_Month <-mydata2[which(mydata2$TSLW<=90*3 & mydata2$TSLW>90*2),]
TSLW_12_Month <-mydata2[which(mydata2$TSLW<=90*4 & mydata2$TSLW>90*3),]
TSLW_1_year<-mydata2[which(mydata2$TSLW>90*4),]

TSLW_1_day=TSLW_1_day[,-c(1:2)]
TSLW_1_Month=TSLW_1_Month[,-c(1:2)]
TSLW_3_Month=TSLW_3_Month[,-c(1,2)]
TSLW_6_Month=TSLW_6_Month[,-c(1,2)]
TSLW_9_Month=TSLW_9_Month[,-c(1,2)]
TSLW_12_Month=TSLW_12_Month[,-c(1,2)]
TSLW_1_year=TSLW_1_year[,-c(1,2)]

TSLW_1_day=TSLW_1_day[rowSums(TSLW_1_day!= 0) > 0,]	# remove columns and or rows with no records
TSLW_1_day=TSLW_1_day[,colSums(TSLW_1_day!= 0) > 0]
TSLW_1_Month=TSLW_1_Month[rowSums(TSLW_1_Month!= 0) > 0,]	# remove columns and or rows with no records
TSLW_1_Month=TSLW_1_Month[,colSums(TSLW_1_Month!= 0) > 0]
TSLW_3_Month=TSLW_3_Month[rowSums(TSLW_3_Month!= 0) > 0,]	# remove columns and or rows with no records
TSLW_3_Month=TSLW_3_Month[,colSums(TSLW_3_Month!= 0) > 0]
TSLW_6_Month=TSLW_6_Month[rowSums(TSLW_6_Month!= 0) > 0,]	# remove columns and or rows with no records
TSLW_6_Month=TSLW_6_Month[,colSums(TSLW_6_Month!= 0) > 0]
TSLW_9_Month=TSLW_9_Month[rowSums(TSLW_9_Month!= 0) > 0,]	# remove columns and or rows with no records
TSLW_9_Month=TSLW_9_Month[,colSums(TSLW_9_Month!= 0) > 0]
TSLW_12_Month=TSLW_12_Month[rowSums(TSLW_12_Month!= 0) > 0,]	# remove columns and or rows with no records
TSLW_12_Month=TSLW_12_Month[,colSums(TSLW_12_Month!= 0) > 0]
TSLW_1_year=TSLW_1_year[rowSums(TSLW_1_year!= 0) > 0,]	# remove columns and or rows with no records
TSLW_1_year=TSLW_1_year[,colSums(TSLW_1_year!= 0) > 0]

TSLW_1_day=as.matrix(TSLW_1_day)
HLWL.PA.1.day<-decostand(TSLW_1_day,method="pa")
HLWL.PA.core.1.day <- betapart.core(HLWL.PA.1.day)
HLWL.PA.multi.1.day <- beta.multi(HLWL.PA.core.1.day)
HLWL.samp.1.day <- beta.sample(HLWL.PA.core.1.day, sites=20, samples=1000)
HLWL.dist.n.1.day <- HLWL.samp.1.day$sampled.values


TSLW_1_Month=as.matrix(TSLW_1_Month)
HLWL.PA.1<-decostand(TSLW_1_Month,method="pa")
HLWL.PA.core.1 <- betapart.core(HLWL.PA.1)
HLWL.PA.multi.1 <- beta.multi(HLWL.PA.core.1)
HLWL.samp.1 <- beta.sample(HLWL.PA.core.1, sites=20, samples=1000)
HLWL.dist.n.1 <- HLWL.samp.1$sampled.values

TSLW_3_Month=as.matrix(TSLW_3_Month)
HLWL.PA.3<-decostand(TSLW_3_Month,method="pa")
HLWL.PA.core.3 <- betapart.core(HLWL.PA.3)
HLWL.PA.multi.3 <- beta.multi(HLWL.PA.core.3)
HLWL.samp.3 <- beta.sample(HLWL.PA.core.3, sites=20, samples=1000)
HLWL.dist.n.3 <- HLWL.samp.3$sampled.values

TSLW_6_Month=as.matrix(TSLW_6_Month)
HLWL.PA.6<-decostand(TSLW_6_Month,method="pa")
HLWL.PA.core.6 <- betapart.core(HLWL.PA.6)
HLWL.PA.multi.6 <- beta.multi(HLWL.PA.core.6)
HLWL.samp.6 <- beta.sample(HLWL.PA.core.6, sites=20, samples=1000)
HLWL.dist.n.6 <- HLWL.samp.6$sampled.values

TSLW_9_Month=as.matrix(TSLW_9_Month)
HLWL.PA.9<-decostand(TSLW_9_Month,method="pa")
HLWL.PA.core.9 <- betapart.core(HLWL.PA.9)
HLWL.PA.multi.9 <- beta.multi(HLWL.PA.core.9)
HLWL.samp.9 <- beta.sample(HLWL.PA.core.9, sites=20, samples=1000)
HLWL.dist.n.9 <- HLWL.samp.9$sampled.values

TSLW_12_Month=as.matrix(TSLW_12_Month)
HLWL.PA.12<-decostand(TSLW_12_Month,method="pa")
HLWL.PA.core.12 <- betapart.core(HLWL.PA.12)
HLWL.PA.multi.12 <- beta.multi(HLWL.PA.core.12)
HLWL.samp.12 <- beta.sample(HLWL.PA.core.12, sites=20, samples=1000)
HLWL.dist.n.12 <- HLWL.samp.12$sampled.values


TSLW_1_year=as.matrix(TSLW_1_year)
HLWL.PA.1yr<-decostand(TSLW_1_year,method="pa")
HLWL.PA.core.1yr <- betapart.core(HLWL.PA.1yr)
HLWL.PA.multi.1yr <- beta.multi(HLWL.PA.core.1yr)
HLWL.samp.1yr <- beta.sample(HLWL.PA.core.1yr, sites=20, samples=1000)
HLWL.dist.n.1yr <- HLWL.samp.1yr$sampled.values
# colours

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 7
cols = gg_color_hue(n)

png(paste(image.dir,"Hattah Lakes TSLW betapart May 2019.png",sep=''),width=24, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(1,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

plot(density(HLWL.dist.n.1.day$beta.SOR), xlim=c(0,1), ylim=c(0, 100), xlab="Beta diversity", main="", lwd=3) 


lines(density(HLWL.dist.n.1.day$beta.SOR), col=cols[1], lwd=3)
lines(density(HLWL.dist.n.1.day$beta.SNE), col=cols[1], lty=1, lwd=2)
lines(density(HLWL.dist.n.1.day$beta.SIM), col=cols[1], lty=2, lwd=2)

lines(density(HLWL.dist.n.1$beta.SOR), col=cols[2], lwd=3)
lines(density(HLWL.dist.n.1$beta.SNE), col=cols[2], lty=1, lwd=2)
lines(density(HLWL.dist.n.1$beta.SIM), col=cols[2], lty=2, lwd=2)

lines(density(HLWL.dist.n.3$beta.SOR),col=cols[3], lwd=3)
lines(density(HLWL.dist.n.3$beta.SNE),col=cols[3],lty=1, lwd=2)
lines(density(HLWL.dist.n.3$beta.SIM),col=cols[3],lty=2, lwd=2)

lines(density(HLWL.dist.n.6$beta.SOR),col=cols[4], lwd=3)
lines(density(HLWL.dist.n.6$beta.SNE),col=cols[4], lty=1, lwd=2)
lines(density(HLWL.dist.n.6$beta.SIM),col=cols[4], lty=2, lwd=2)

lines(density(HLWL.dist.n.9$beta.SOR),col=cols[5], lwd=3)
lines(density(HLWL.dist.n.9$beta.SNE),col=cols[5], lty=1, lwd=2)
lines(density(HLWL.dist.n.9$beta.SIM),col=cols[5], lty=2, lwd=2)

lines(density(HLWL.dist.n.12$beta.SOR),col=cols[6], lwd=3)
lines(density(HLWL.dist.n.12$beta.SNE),col=cols[6], lty=1, lwd=2)
lines(density(HLWL.dist.n.12$beta.SIM),col=cols[6], lty=2, lwd=2)

lines(density(HLWL.dist.n.1yr$beta.SOR),col=cols[7], lwd=3)
lines(density(HLWL.dist.n.1yr$beta.SNE),col=cols[7], lty=1, lwd=2)
lines(density(HLWL.dist.n.1yr$beta.SIM),col=cols[7], lty=2, lwd=2)


legend(0, 100, legend=c("Inundated","1 Month", "3 Months","6 Months","9 Months", "12 Months",">1 year"),col=cols, lty=1, lwd=2,cex=1)

dev.off()
