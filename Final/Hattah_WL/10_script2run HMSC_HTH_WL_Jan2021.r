####=============================================================================================================================
####Code to run HMSC on Hattah data
# Import data

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
mydata=read.csv("Spp_site_year_transect matrix HTH_WL_May 2019.csv",row.names=1) # save data out

date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
env.data<-read.csv("Hattah Lakes wetlands transect based env data January 2021.csv") # 

# remove duplicated rows

env.data2=env.data[!duplicated(env.data$Unique_site_year_season),]

# Tidy up environmental data and remove rows where the metric can not be calculated

envdata=env.data2
envdata=merge(envdata,mydata,by.x="Unique_site_year_season",by.y="row.names")


####=============================================================================================================================
# set up and run Hmsc model
library(Hmsc)


Y=envdata[,41:250] # subset to species only			
rownames(Y) =envdata$Unique_site_year_season																			# remove grouping column (1) as now row names
YY=Y
Y[Y>0]<-1 																										# turn data to presence/absence
sp_occur<-apply(Y,2,sum)																							# Sum species occurrences 
mydat_common =YY[,sp_occur>4] 		
mypreds=envdata[,1:40]
Alldat=merge(mypreds, mydat_common, by.x="Unique_site_year_season", by="row.names")
AlldatV1 = Alldat[rowSums(Alldat[,41:155])!=0, ] # 
Y=AlldatV1[,41:155]

AlldatV1$d90 <-scale(AlldatV1$d90)
AlldatV1$d365 <-scale(AlldatV1$d365)
AlldatV1$TSLW_LT<-scale(AlldatV1$TSLW_LT)
AlldatV1$MeanTemp90 <-scale(AlldatV1$MeanTemp90)
AlldatV1$MaxTemp90 <-scale(AlldatV1$MaxTemp90)
AlldatV1$MeanTemp365 <-scale(AlldatV1$MeanTemp365)
AlldatV1$MaxTemp365 <-scale(AlldatV1$MaxTemp365)
AlldatV1$Wetland <-as.factor(AlldatV1$Wetland)
AlldatV1$VEG_CLASS <-as.factor(AlldatV1$VEG_CLASS)
AlldatV1$Inundated  <-as.factor(AlldatV1$Inundated )


Y=as.matrix(Y) 
XData=data.frame(AlldatV1[,1:40])
xy=as.matrix(cbind(AlldatV1$Easting, AlldatV1$Northing))
XData=data.frame(x1=XData$TSLW_LT)
Xmatrix=as.matrix(XData)
studyDesign=data.frame(Plot=droplevels(as.factor(AlldatV1$Wetland)),Site=droplevels(as.factor(AlldatV1$Site.ID)))
rL1 =HmscRandomLevel(units = studyDesign$Plot) # set a random effect for waterhole
rL2 =HmscRandomLevel(units = studyDesign$Site) # set a random effect for site
rL3 =HmscRandomLevel(sData=xy)# add spatial random effect

XData=data.frame(AlldatV1[,1:40])
xy=as.matrix(cbind(AlldatV1$Easting, AlldatV1$Northing))
XData=data.frame(x1=XData$TSLW_LT)
XFormula=~poly(x1, degree=2, raw = TRUE)
m_TSLW=Hmsc(Y=Y, XData=XData, XFormula=XFormula, distr=c("probit"),studyDesign=studyDesign, ranLevels=list(route=rL3, Plot=rL1, Site=rL2))

XData=data.frame(AlldatV1[,1:40])
xy=as.matrix(cbind(AlldatV1$Easting, AlldatV1$Northing))
XData=data.frame(x1=XData$TSLW_LT, x2=XData$d3Mon_wet, x3=XData$d3Mon_meandepth)

XFormula=~poly(x1, degree=2, raw = TRUE)+x2+x3
m_TSLW_d3Mon=Hmsc(Y=Y, XData=XData, XFormula=XFormula, distr=c("probit"),studyDesign=studyDesign, ranLevels=list(route=rL3, Plot=rL1, Site=rL2))

XData=data.frame(AlldatV1[,1:40])
xy=as.matrix(cbind(AlldatV1$Easting, AlldatV1$Northing))
XData=data.frame(x1=XData$TSLW_LT, x2=XData$d3Mon_wet, x3=XData$d3Mon_meandepth, x4=XData$d1yrs_wet, x5=XData$d1yrs_meandepth)

XFormula=~poly(x1, degree=2, raw = TRUE)+x2+x3+x4+x5
m_TSLW_d1yr=Hmsc(Y=Y, XData=XData, XFormula=XFormula, distr=c("probit"),studyDesign=studyDesign, ranLevels=list(route=rL3, Plot=rL1, Site=rL2))

XData=data.frame(AlldatV1[,1:40])
xy=as.matrix(cbind(AlldatV1$Easting, AlldatV1$Northing))
XData=data.frame(x1=XData$TSLW_LT, x2=XData$d3Mon_wet, x3=XData$d3Mon_meandepth, x4=XData$d1yrs_wet, x5=XData$d1yrs_meandepth,x6=XData$d3yrs_wet,x7=XData$d3yrs_wet)

XFormula=~poly(x1, degree=2, raw = TRUE)+x2+x3+x4+x5+x6+x7
m_TSLW_d3yr=Hmsc(Y=Y, XData=XData, XFormula=XFormula, distr=c("probit"),studyDesign=studyDesign, ranLevels=list(route=rL3, Plot=rL1, Site=rL2))

XData=data.frame(AlldatV1[,1:40])
xy=as.matrix(cbind(AlldatV1$Easting, AlldatV1$Northing))
XData=data.frame(x1=XData$TSLW_LT, x2=XData$Max.CTF.30yrs, x3=XData$n.events.30yrs, x4=XData$FF.30, x5=XData$P.CTF.30yrs)

XFormula=~poly(x1, degree=2, raw = TRUE)+x2+x3+x4+x5
m_TSLW_d30yr=Hmsc(Y=Y, XData=XData, XFormula=XFormula, distr=c("probit"),studyDesign=studyDesign, ranLevels=list(route=rL3, Plot=rL1, Site=rL2))

XData=data.frame(AlldatV1[,1:40])
xy=as.matrix(cbind(AlldatV1$Easting, AlldatV1$Northing))
XData=data.frame(x1=XData$TSLW_LT, x2=XData$Max.CTF.10yrs, x3=XData$n.events.10yrs, x4=XData$FF.10,x5=XData$P.CTF.10yrs)

XFormula=~poly(x1, degree=2, raw = TRUE)+x2+x3+x4+x5
m_TSLW_d10yr=Hmsc(Y=Y, XData=XData, XFormula=XFormula, distr=c("probit"),studyDesign=studyDesign, ranLevels=list(route=rL3, Plot=rL1, Site=rL2))


nChains = 2
test.run = FALSE
if(test.run){#with this option, the vignette runs fast but results are not reliable
thin = 1
samples = 10

transient = 5
verbose = 5}else{#with this option, the vignette evaluates slow but it reproduces the results of the#.pdf version
thin = 35
samples = 1000
transient = 15000
verbose = 15000
}

# run without climate
models=list(m_TSLW, m_TSLW_d3Mon,m_TSLW_d1yr, m_TSLW_d3yr, m_TSLW_d10yr, m_TSLW_d30yr)

for(i in 1:6){

models=sampleMcmc(models[[i]], thin=thin, samples=samples, transient=transient, nChains=nChains, verbose=verbose, initPar="fixed effects")
}


m =sampleMcmc(m, thin = thin, samples = samples, transient = transient,nChains = nChains, verbose = verbose)


# explore convergence
mpost=convertToCodaObject(m)
effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta[,1:10])
postBeta = getPostEstimate(m, parName = "Beta")
plot(mpost$Beta[,1:2])
plotBeta(m,
         post = postBeta, 
         plotTree = F,
         spNamesNumbers = c(T,T))
		 
#explore fit	 
		
preds=computePredictedValues(m)		
Mfit=evaluateModelFit(hM=m, predY=preds)
Mfit$R2

save(m, file = "Probit Hmsc model 20200608.RData")

load("C:/Users/jc246980/Documents/Current projects/TFTA/TFTA invert analysis/Lognormal poisson Hmsc model 20200316.RData")

####=============================================================================================================================
# explore mcmc convergence etc
mpost=convertToCodaObject(m)

# in examining diagnotistics for residual species associations 48 species means that the matrix is huge so we select a random subsample of 100 sleected species pairs to explore

sppairs=matrix(sample(x=1:64^2,size=100))
tmp = mpost$Omega[[1]]

for(chain in 1:length(tmp)){
tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega =effectiveSize(tmp)
psrf.omega =gelman.diag(tmp, multivariate=FALSE)$psrfpng(paste(image.dir,"TFTA_HMSC LN poisson diagnostics 20200319.png",sep=''),width=20, height=20, units='cm', res=400, pointsize=8, bg='white')
par(mfrow=c(2,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf,main="psfr(beta)")
hist(effectiveSize(mpost$Omega[[1]]), main="ess(omega)")
hist(psrf.omega,main="psfr(omega) subsample 100 pairs")



dev.off()

####=============================================================================================================================
# cross validation over 51 separate plots (leave one 'plot' out for each run) - not sure we want to do this for full model as it might take a year to run!!

partition=createPartition(m,nfolds=17,column="Plot")
preds=computePredictedValues(m, partition=partition)
MF=evaluateModelFit(hM=m, predY=preds)

####=============================================================================================================================
# evalulate model fit, explanatory power

preds=computePredictedValues(m, expected=FALSE)
MF=evaluateModelFit(hM=m, predY=preds)

myeval <-data.frame(MF$RMSE, MF$SR2, MF$O.AUC, MF$O.TjurR2, MF$O.RMSE, MF$C.SR2, MF$C.RMSE)
write.csv(myeval , file = "TFTA_HMSC_model evaluations 20200319.csv") #

####=============================================================================================================================
# get posterior distribution - means and quantiles (credibility interval) for each species

postBeta=getPostEstimate(m, parName="Beta",q=c(0.025,0.5,0.975))

png(paste(image.dir,"TFTA_HMSC LN means BetaPlot 202010528.png",sep=''),width=18, height=30, units='cm', res=400, pointsize=8, bg='white')
    par(mar=c(6,6,4,2),cex=0.7) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns	

plotBeta(m,post=postBeta, para="Support", supportLevel=0.95)

dev.off()

Support <- postBeta$support
SupportNeg <- postBeta$supportNeg

SupportALL <-rbind(Support, SupportNeg)

write.csv(SupportALL , file = "TFTA_HMSC_SupportFiles 20200319.csv")

mybetas=as.matrix(summary(mpost$Beta))
write.csv(mybetas, file = "TFTA_HMSC_beta estimates 20200319.csv")

####=============================================================================================================================
# visualising the residuals between species
library(corrplot)

png(paste(image.dir,"TFTA_HMSC residual correlation between sp.png",sep=''),width=20, height=20, units='cm', res=400, pointsize=8, bg='white')
    par(mar=c(6,6,4,2),cex=0.7) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns	
	
OmegaCor=computeAssociations(m)
supportLevel=0.95
toPlot=((OmegaCor[[1]]$support>supportLevel)+
(OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
corrplot(toPlot,method="color",
col=colorRampPalette(c("blue", "white", "red"))(200),
title=paste("random effect level:", m$rLNames[1]), mar=c(0,0,1,0))

dev.off()

write.csv(OmegaCor , file = "TFTA_HMSC_Residual sp correlations.csv") #

####=============================================================================================================================
# Variance partitioning

head(m$X)

XData=data.frame(x1=XData$Habitat, x2=XData$Rain6M,x3=XData$Rain18M,x4=XData$TSS,x5=XData$DIN,x6=XData$TN, 
x7=XData$FRP, x8=XData$TP,x9=XData$Cond,x10=XData$TempMax,x11=XData$TSLF50, x12=XData$TSLF25)

VP =computeVariancePartitioning(m, group =c(1,2,2,3,3,3,3,3,3,4,3,3), groupnames =c("Habitat", "Rain","WQ","Temp"))


png(paste(image.dir,"TFTA_HMSC_variance_partioning 202020319.png",sep=''),width=30, height=15, units='cm', res=400, pointsize=8, bg='white')
    par(mar=c(6,6,4,2),cex=0.7) #defines plot parameters mfrow is number of images per sheet in this case its 5 rows by two columns	
	
plotVariancePartitioning(m, VP = VP, las=2)

dev.off()

write.csv(VP$vals, file = "TFTA_HMSC_variance partioning 20200319.csv")

####=============================================================================================================================
# Exploring responses to gradients of continuous predictors

Gradient=constructGradient(m,focalVariable="x9",ngrid=20)
predY=predict(m, XData=Gradient$XDataNew, studyDesign=Gradient$studyDesignNew, ranLevels=Gradient$rLNew, expected=TRUE)
par(mfrow=c(3,3))
plotGradient(m, Gradient, pred=predY, measure="S", las=1,showData=TRUE, main="Species richness")
plotGradient(m, Gradient, pred=predY, measure="Y", index=3, las=1,showData=TRUE, main="Focal species")

####=============================================================================================================================