### Script to determine beta diversity calculations for each wetland
### Cassie James
### 8 June 2018

library(betapart)
library(vegan)

#############################################################################################
### Code to determine beta diversity for all wetland complexes together
image.dir="C:/Users/jc246980/Documents/MD Vegetation/Plots/"

data.dir="C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.HLWL=read.csv("Spp_site matrix HTH_WL_July 2018.csv",row.names=1) # load data - object name is 'tada' ... :)
Output.HLWL =data.matrix.HLWL
Output.HLWL=Output.HLWL[rowSums(Output.HLWL!= 0) > 0,]	# remove columns and or rows with no records
Output.HLWL=Output.HLWL[,colSums(Output.HLWL!= 0) > 0]

data.dir="C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd(data.dir)
data.matrix.LMW=read.csv("Spp_site matrix LMW_WL_July 2018.csv", row.names=1) # load data - object name is 'Output' ... :)
Output.LMW =data.matrix.LMW
Output.LMW=Output.LMW[rowSums(Output.LMW!= 0) > 0,]	# remove columns and or rows with no records
Output.LMW=Output.LMW[,colSums(Output.LMW!= 0) > 0]

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
data.matrix.Gun=read.csv("Spp_site_matrix Gunbower_WL_July 2018.csv",row.names=1) # load data - object name is 'tada' ... :)
Output.Gun =data.matrix.Gun
Output.Gun=Output.Gun[rowSums(Output.Gun!= 0) > 0,]	# remove columns and or rows with no records
Output.Gun=Output.Gun[,colSums(Output.Gun!= 0) > 0]
Output.Gun=Output.Gun[,-grep("NA",colnames(Output.Gun))] # remove 'NA' column

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/KP_data_csvs/"; setwd (data.dir) 
data.matrix.KP=read.csv("Spp_site_matrix KP_WL_July 2018.csv",row.names=1) # load data - object name is 'tada' ... :)
Output.KP =data.matrix.KP
Output.KP=Output.KP[,-grep("NA",colnames(Output.KP))] # remove 'NA' column
Output.KP=Output.KP[rowSums(Output.KP!= 0) > 0,]	# remove columns and or rows with no records
Output.KP=Output.KP[,colSums(Output.KP!= 0) > 0]

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Barmah_data_csvs/"; setwd (data.dir) 
data.matrix.Barm=read.csv("Spp_site_matrix Barmah_WL_July 2018.csv",row.names=1) # load data - object name is 'tada' ... :)
Output.Barm =data.matrix.Barm
Output.Barm=Output.Barm[rowSums(Output.Barm!= 0) > 0,]	# remove columns and or rows with no records
Output.Barm=Output.Barm[,colSums(Output.Barm!= 0) >0]

data.dir="C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd(data.dir)
data.matrix.Chow=read.csv("Spp_site_matrix Chowilla_WL_July 2018.csv",row.names=1) # load data - object name is 'tada' ... :)
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

png(paste(image.dir,"Wetlands together betapart.png",sep=''),width=24, height=12, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(1,1),cex=1,oma=c(2,0,1,0.5)) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns

plot(density(HLWL.dist.n$beta.SOR), xlim=c(0,0.8), ylim=c(0, 70), xlab="Beta diversity", main="", lwd=3) 

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

legend(0, 70, legend=c("Barmah", "Chowilla","Gunbower","Hattah", "KP", "LMW"),col=cols, lty=1, lwd=2,cex=1)

dev.off()

##############################################################################################################
#### Code to look at Hattah lakes wetlands temporal beta diversity and spatial beta diversity

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
