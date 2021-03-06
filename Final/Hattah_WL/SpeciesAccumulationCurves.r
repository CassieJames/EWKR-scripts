##########################################################################################################################################
#### Script for species accumulation curves for Hattah wetlands
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"
data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir) 
Output=read.csv("Spp_site_year_transect matrix HTH_WL_May 2019 ABUND.csv", row.names=1) # import 


rownames (Output) <- gsub("LHT", "LHAT", rownames(Output)) 
rownames (Output) <- gsub("CHT", "CCS", rownames(Output)) 

wetlist=c("BIT","BLT", "BOT", "BRT", "CCS", "HT","LHAT", "KT", "MOT", "NCT", "NN", "YT")

png(paste(image.dir,"Hattah Wetlands accumulation_curves_May 2019.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(4,3),cex=1,oma=c(2,0,1,0.5))
		
for(w in wetlist) { # 
tdata=Output[grep(w,rownames(Output)),]
tdata[tdata>0] <-1

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

set.seed(101)
accum <- specaccum(tdata, method="exact")
slopes <- with(accum,diff(richness)/diff(sites))
plat <-which(slopes<0.1)[1]
plateau<-round(accum$richness[plat],0)
Exact<-accum$richness[length(accum$richness)]

sp1 <- specaccum(tdata, method = "collector")
MM <- fitspecaccum(accum,  "michaelis-menten")
assym <- fitspecaccum(accum,  "asymp")
chao <-round(specpool(tdata)$chao,0)
boot <-round(specpool(tdata)$boot,0)

z <- betadiver(tdata, "z")

plot(accum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlab=w, ylab="Cumulative species number" )
plot(sp1, add=TRUE,col="red")

legend('bottomright', legend = c(paste("Chao =", chao),paste("Exact =",Exact)))

}
dev.off()

##########################################################################################################################################
#### Script for species accumulation curves for LMW wetlands

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 
Output=read.csv("Spp_site_year_transect matrix LMW_May 2019 ABUND.csv", row.names=1) # import data

wetlist=c("BB","CR","LP","UL","MUH","BI","UMWC","W33","SCB","MLH","WL","WW")

png(paste(image.dir,"LMW Wetlands accumulation_curves_May 2019.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(4,3),cex=1,oma=c(2,0,1,0.5))
		
for(w in wetlist) { # 
tdata=Output[grep(w,rownames(Output)),]
tdata[tdata>0] <-1

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

set.seed(101)
accum <- specaccum(tdata, method="exact")
slopes <- with(accum,diff(richness)/diff(sites))
plat <-which(slopes<0.1)[1]
plateau<-round(accum$richness[plat],0)
Exact<-accum$richness[length(accum$richness)]

sp1 <- specaccum(tdata, method = "collector")
MM <- fitspecaccum(accum,  "michaelis-menten")
assym <- fitspecaccum(accum,  "asymp")
chao <-round(specpool(tdata)$chao,0)
boot <-round(specpool(tdata)$boot,0)

z <- betadiver(tdata, "z")

plot(accum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlab=w, ylab="Cumulative species number" )
plot(sp1, add=TRUE,col="red")

legend('bottomright', legend = c(paste("Chao =", chao),paste("Exact =",Exact)))

}
dev.off()



##########################################################################################################################################
#### Script for species accumulation curves for Gunbower wetlands

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
Output=read.csv("Spp_site_transect_year_season matrix Gunbower_WL_May 2019 ABUND.csv", row.names=1) # import data

wetlist=c("LL","GS","LG","IP","RL","BLS","FB","CS","LR","COS")


png(paste(image.dir,"Gunbower Wetlands accumulation_curves_May 2019.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(4,3),cex=1,oma=c(2,0,1,0.5))
		
for(w in wetlist) { # 
tdata=Output[grep(w,rownames(Output)),]
tdata[tdata>0] <-1

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

set.seed(101)
accum <- specaccum(tdata, method="exact")
slopes <- with(accum,diff(richness)/diff(sites))
plat <-which(slopes<0.1)[1]
plateau<-round(accum$richness[plat],0)
Exact<-accum$richness[length(accum$richness)]

sp1 <- specaccum(tdata, method = "collector")
MM <- fitspecaccum(accum,  "michaelis-menten")
assym <- fitspecaccum(accum,  "asymp")
chao <-round(specpool(tdata)$chao,0)
boot <-round(specpool(tdata)$boot,0)

z <- betadiver(tdata, "z")

plot(accum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlab=w, ylab="Cumulative species number" )
plot(sp1, add=TRUE,col="red")

legend('bottomright', legend = c(paste("Chao =", chao),paste("Exact =",Exact)))

}
dev.off()

##########################################################################################################################################
#### Script for species accumulation curves for KP wetlands

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/KP_data_csvs/"; setwd (data.dir) 
Output=read.csv("Spp_site_year_transect_matrix_KP_WL_May 2019.csv", row.names=1) # import data

wetlist=c("PR1","SL","WH","BW","CLT","TL","PS","PAW","BC","PRW","PB","PJW") # need to make sure that the wetland names are not subsets of each ot

png(paste(image.dir,"KP Wetlands accumulation_curves_method May 2019.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(4,3),cex=1,oma=c(2,0,1,0.5))
		
for(w in wetlist) { # 
tdata=Output[grep(w,rownames(Output)),]
tdata[tdata>0] <-1

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

set.seed(101)
accum <- specaccum(tdata, method="exact")
slopes <- with(accum,diff(richness)/diff(sites))
plat <-which(slopes<0.1)[1]
plateau<-round(accum$richness[plat],0)
Exact<-accum$richness[length(accum$richness)]

sp1 <- specaccum(tdata, method = "collector")
MM <- fitspecaccum(accum,  "michaelis-menten")
assym <- fitspecaccum(accum,  "asymp")
chao <-round(specpool(tdata)$chao,0)
boot <-round(specpool(tdata)$boot,0)

z <- betadiver(tdata, "z")

plot(accum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlab=w, ylab="Cumulative species number" )
plot(sp1, add=TRUE,col="red")

legend('bottomright', legend = c(paste("Chao =", chao),paste("Exact =",Exact)))

}
dev.off()


wetlist=c("BL","PLL") 

png(paste(image.dir,"KP Wetlands accumulation_curves_May 2019 Part 2.png",sep=''),width=20, height=25, units='cm', res=500, pointsize=10, bg='white')
        par(mar=c(4,4,1,1),mfrow=c(4,3),cex=1,oma=c(2,0,1,0.5))
		
for(w in wetlist) { # 
tdata=Output[grep(w,rownames(Output)),]
tdata[tdata>0] <-1

tdata=tdata[,colSums(tdata!= 0) > 0]	
tdata=tdata[rowSums(tdata!= 0) > 0,]	

set.seed(101)
accum <- specaccum(tdata, method="exact")
slopes <- with(accum,diff(richness)/diff(sites))
plat <-which(slopes<0.1)[1]
plateau<-round(accum$richness[plat],0)
Exact<-accum$richness[length(accum$richness)]

sp1 <- specaccum(tdata, method = "collector")
MM <- fitspecaccum(accum,  "michaelis-menten")
assym <- fitspecaccum(accum,  "asymp")
chao <-round(specpool(tdata)$chao,0)
boot <-round(specpool(tdata)$boot,0)

z <- betadiver(tdata, "z")

plot(accum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlab=w, ylab="Cumulative species number" )
plot(sp1, add=TRUE,col="red")

legend('bottomright', legend = c(paste("Chao =", chao),paste("Exact =",Exact)))

}
dev.off()