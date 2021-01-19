####=============================================================================================================================
####Code to run HMSC on Hattah data
# Import data

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
mydata=read.csv("Spp_site_year_transect matrix HTH_WL_May 2019.csv",row.names=1) # save data out

date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
env.data<-read.csv("Hattah Lakes wetlands transect based env data.csv") # 

# remove duplicated rows

env.data2=env.data[!duplicated(env.data$Unique_site_year_season),]

# Tidy up environmental data and remove rows where the metric can not be calculated

envdata=env.data2
envdata$FF=envdata$d30yearsFF+envdata$Freq_ALL
envdata=envdata[!is.na(envdata$TSLW),] # remove sites where extent of TSLW is beyond hydrological record
envdata$TSLW[envdata$Innundated==TRUE]<-1 # if site is recorded as inundated TRUE change TSLW to 1 day
envdata$d3Mon_wet[is.na(envdata$d3Mon_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d1yrs_wet[is.na(envdata$d1yrs_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d3yrs_wet[is.na(envdata$d3yrs_wet)]<-0 # if site is never wet during time period replace with 0
envdata$d3Mon_meandepth[is.na(envdata$d3Mon_meandepth)]<-0 # if site is never wet during time period replace with 0
envdata$d1yrs_meandepth[is.na(envdata$d1yrs_meandepth)]<-0 # if site is never wet during time period replace with 0
envdata$d3yrs_meandepth[is.na(envdata$d3yrs_meandepth)]<-0 # if site is never wet during time period replace with 0

envdata=envdata[,c("Unique_site_year", "Site.ID.x","Unique_site_year_season","Date.of.collection.x","Elevation", "Easting", "Northing", "Season","Wetland","d90", "d365","TSLW","TSLW_LT","d3Mon_wet",
"d3Mon_meandepth","d3Mon_drylength", "d1yrs_wet", "d1yrs_meandepth","d1yrs_drylength", "d3yrs_wet", "d3yrs_meandepth","d3yrs_drylength","Freq_d1", "Freq_d3", "Freq_d5", "Freq_d10", "Freq_ALL", "WaterYr", 
"Inundated", "d5yearsFF", "d10yearsFF", "d20yearsFF", "d30yearsFF", "MeanTemp90", "MaxTemp90","MinTemp90", "MeanTemp365", "MaxTemp365", "MinTemp365",
"VEG_CLASS", "Max.CTF.30yrs", "P.CTF.30yrs","Max.CTF.20yrs", "P.CTF.20yrs","Max.CTF.10yrs", "P.CTF.10yrs","Max.CTF.5yrs", "P.CTF.5yrs")]


envdata=merge(envdata,mydata,by.x="Unique_site_year_season",by.y="row.names")


####=============================================================================================================================
# set up and run Hmsc model
library(Hmsc)


Y=envdata[,49:258] # subset to species only			
rownames(Y) =envdata$Unique_site_year_season																			# remove grouping column (1) as now row names
YY=Y
Y[Y>0]<-1 																										# turn data to presence/absence
sp_occur<-apply(Y,2,sum)																							# Sum species occurrences 
mydat_common =YY[,sp_occur>9] 		

mypreds=envdata[,1:48]

Alldat=merge(mypreds, mydat_common, by.x="Code", by="row.names")

AlldatV1 = Alldat[rowSums(Alldat[,24:89])!=0, ] # 


Y=AlldatV1[,27:89]

AlldatV1$Rain6M <-scale(AlldatV1$Rain6M)
AlldatV1$Rain18M <-scale(AlldatV1$Rain18M)
AlldatV1$TSS <-scale(AlldatV1$TSS)
AlldatV1$DIN <-scale(AlldatV1$DIN)
AlldatV1$TN <-scale(AlldatV1$TN)
AlldatV1$FRP <-scale(AlldatV1$FRP)
AlldatV1$TP <-scale(AlldatV1$TP)
AlldatV1$Cond <-scale(AlldatV1$Cond)
AlldatV1$Alk <-scale(AlldatV1$Alk)
AlldatV1$LabpH <-scale(AlldatV1$LabpH)
AlldatV1$TChl <-scale(AlldatV1$TChl)
AlldatV1$Hard <-scale(AlldatV1$Hard)
AlldatV1$TempMax <-scale(AlldatV1$TempMax)
AlldatV1$TSLF50 <-scale(AlldatV1$TSLF50)
AlldatV1$TSLF25 <-scale(AlldatV1$TSLF25)

Y=as.matrix(Y) # species matrix
XData=data.frame(AlldatV1[,1:26])


XData=data.frame(x1=XData$Habitat, x2=XData$Rain6M,x3=XData$Rain18M,x4=XData$TSS,x5=XData$DIN,x6=XData$TN, 
x7=XData$FRP, x8=XData$TP,x9=XData$Cond,x10=XData$TempMax,x11=XData$TSLF50, x12=XData$TSLF25, x13=XData$Hard, x14=XData$Alk, x15=XData$LabpH, x16=XData$TChl)

Xmatrix=as.matrix(XData)

studyDesign=data.frame(Sample=droplevels(as.factor(AlldatV1$Code)), Plot=droplevels(as.factor(AlldatV1$Site)))
rL1 =HmscRandomLevel(units = studyDesign$Plot) # set a random effect for waterhole
rL2 =HmscRandomLevel(units = studyDesign$Sample) # set a random effect at the level of sampling unit


m=Hmsc(Y=Y, XData=XData, XFormula=~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16, distr=c("lognormal poisson"),studyDesign=studyDesign, ranLevels=list(Sample=rL2, Plot=rL1))

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
verbose = 50*thin
}

m =sampleMcmc(m, thin = thin, samples = samples, transient = transient,nChains = nChains, verbose = verbose)
save(m, file = "Lognormal poisson Hmsc model 20200608.RData")

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