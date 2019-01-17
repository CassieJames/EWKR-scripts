##############################################################################
#### Model richness and abundance of wetland and dryland species 
library(corrplot)
library(mboost)


data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Hattah wetlands response by metrics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"

# sort out environmental data - centre and log
envdata=data.matrix.env.data
envdata$TSLW[envdata$Inundated==TRUE] <- 1 # If site is currently inundated make TSLW 1 day

envdata$Inundated = as.factor(envdata$Inundated )
# highly skewed variables were log10 transformed before analysis
envdata$d30=as.numeric(scale(log10(envdata$d30+1),center=TRUE, scale=FALSE)) # 
envdata$d90=as.numeric(scale(log10(envdata$d90+1),center=TRUE, scale=FALSE)) # added one to avoid return of infinity values due to low mean rainfall over shorter time periods
envdata$d365=as.numeric(scale(log10(envdata$d365+1),center=TRUE, scale=FALSE))
#envdata$TSLW=as.numeric(scale(log10(envdata$TSLW+1),center=TRUE, scale=FALSE)) this is a very strongly skewed variable as there are lots of '1's which indicate its still wte
envdata$MRI_shallow=as.numeric(scale(log10(envdata$MRI_shallow+1),center=TRUE, scale=FALSE)) # Highly skewed - log helps a bit but its heavily weighted to the left
envdata$MRI_medium=as.numeric(scale(log10(envdata$MRI_medium+1),center=TRUE, scale=FALSE))
envdata$MRI_deep=as.numeric(scale(log10(envdata$MRI_deep+1),center=TRUE, scale=FALSE))
envdata$d1yrs_dry=as.numeric(scale((envdata$d1yrs_dry),center=TRUE, scale=FALSE)) # quite strong peaks at low and high end of range
envdata$d1yrs_wet=as.numeric(scale((envdata$d1yrs_wet),center=TRUE, scale=FALSE)) # quite strong peaks at low and high end of range
envdata$d1yrs_shallow=as.numeric(scale(log10(envdata$d1yrs_shallow+1),center=TRUE, scale=FALSE)) #this is a very strongly skewed variable as there are lots of low values
envdata$d1yrs_medium=as.numeric(scale(log10(envdata$d1yrs_medium+1),center=TRUE, scale=FALSE))#this is a very strongly skewed variable as there are lots of low values
envdata$d1yrs_deep=as.numeric(scale(log10(envdata$d1yrs_deep+1),center=TRUE, scale=FALSE))
envdata$d3yrs_dry=as.numeric(scale((envdata$d3yrs_dry),center=TRUE, scale=FALSE))
envdata$d3yrs_shallow=as.numeric(scale(envdata$d3yrs_shallow,center=TRUE, scale=FALSE)) #this is a very strongly skewed variable as there are lots of low values
envdata$d3yrs_medium=as.numeric(scale(envdata$d3yrs_medium,center=TRUE, scale=FALSE))
envdata$d3yrs_deep=as.numeric(scale(envdata$d3yrs_deep,center=TRUE, scale=FALSE))
envdata$Easting=as.numeric(scale(envdata$Easting^2,center=FALSE, scale=TRUE))
envdata$Northing=as.numeric(scale(envdata$Northing^2,center=FALSE, scale=TRUE))
envdata$INT <- rep(1,nrow(envdata)) # provide intercept variable
envdata$MaxTemp30=as.numeric(scale((envdata$MaxTemp30),center=TRUE, scale=FALSE))
envdata$MinTemp30=as.numeric(scale((envdata$MinTemp30),center=TRUE, scale=FALSE))
envdata$MeanTemp30=as.numeric(scale((envdata$MeanTemp30),center=TRUE, scale=FALSE))
envdata$MaxTemp90=as.numeric(scale((envdata$MaxTemp90),center=TRUE, scale=FALSE))
envdata$MinTemp90=as.numeric(scale((envdata$MinTemp90),center=TRUE, scale=FALSE))
envdata$MeanTemp90=as.numeric(scale((envdata$MeanTemp90),center=TRUE, scale=FALSE))
envdata$MaxTemp365=as.numeric(scale((envdata$MaxTemp365),center=TRUE, scale=FALSE))
envdata$MinTemp365=as.numeric(scale((envdata$MinTemp365),center=TRUE, scale=FALSE))
envdata$MeanTemp365=as.numeric(scale((envdata$MeanTemp365),center=TRUE, scale=FALSE))

envdata.short <- envdata[,c("d30","d90","d365", "TSLW", "MRI_shallow", "MRI_medium","MRI_deep","d1yrs_dry","d1yrs_shallow","d1yrs_medium","d1yrs_deep","d3yrs_dry","d3yrs_shallow","d3yrs_medium","d3yrs_deep",
"MaxTemp30","MinTemp30", "MeanTemp30","MaxTemp90","MinTemp90", "MeanTemp90","MaxTemp365","MinTemp365", "MeanTemp365")]

envdata.short.matrix<-as.matrix(as.data.frame(envdata.short)) 
M <- cor(envdata.short.matrix)
corrplot(M,type = "upper", tl.col = "black", tl.srt = 45) 

daten=envdata

daten$Wet_NativesPA <- 0
daten$Wet_NativesPA[daten$Wet_Natives > 0] <-1
daten$Wet_NativesPA=as.factor(daten$Wet_NativesPA)

daten$WetNatRichPA <- 0
daten$WetNatRichPA[daten$WetNatRich > 0] <-1
daten$WetNatRichPA=as.factor(daten$WetNatRichPA)

#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- WetNatRichPA ~ bols(Easting, intercept=FALSE)+bols(Northing, intercept=FALSE)+
bspatial(Easting, Northing, knots = 20, boundary.knots=NULL,center=TRUE, df=1, differences=1)+
bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=1)+
bols(MRI_shallow, intercept=FALSE)+bbs(MRI_shallow, center=TRUE, df=1)+
bols(MRI_medium, intercept=FALSE)+bbs(MRI_medium, center=TRUE, df=1)+
bols(MRI_deep, intercept=FALSE)+bbs(MRI_deep, center=TRUE, df=1)+
bols(d1yrs_dry, intercept=FALSE)+bbs(d1yrs_dry, center=TRUE, df=1)+
bols(d1yrs_shallow, intercept=FALSE)+bbs(d1yrs_shallow, center=TRUE, df=1)+
bols(d1yrs_medium, intercept=FALSE)+bbs(d1yrs_medium, center=TRUE, df=1)+
bols(d1yrs_deep, intercept=FALSE)+bbs(d1yrs_deep, center=TRUE, df=1)+
bols(d3yrs_dry, intercept=FALSE)+bbs(d3yrs_dry, center=TRUE, df=1)+
bols(d3yrs_shallow, intercept=FALSE)+bbs(d3yrs_shallow, center=TRUE, df=1)+
bols(d3yrs_medium, intercept=FALSE)+bbs(d3yrs_medium, center=TRUE, df=1)+
bols(d3yrs_deep, intercept=FALSE)+bbs(d3yrs_deep, center=TRUE, df=1)+
bols(Inundated, intercept=FALSE)


#specify model formula - kept relatively simple of the time being (bols are linear effects, bbs are smoothed effects and bspatial are spatial effects)
formulaB <- WetNatRichPA ~ bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=1)+brandom(Site.ID,df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=1)+
bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=1)+
bols(d3yrs_dry, intercept=FALSE)+bbs(d3yrs_dry, center=TRUE, df=1)+
bols(Inundated, intercept=FALSE)


# Create test and training data sets

n<- nrow(daten) # 
set.seed(806)
indvecL <-sample(1:n,n,replace=FALSE) # create a random set of numbers
datenSmall <- daten[indvecL,][1:ceiling(n/1.7),] # create a subset of data (70%) for training set
test_start=ceiling(n/1.7)+1
testdata <- daten[indvecL,][test_start:length(indvecL),] # create a subset of data (30%) for test set

#data = subset(datenSmall, Wet_Natives > 0),

# Run model on training set
Full.model <-gamboost(formulaB,family = Binomial(), data = datenSmall, control=boost_control(mstop=1000,trace=TRUE)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if intercept =False for Easting, North


formulaB <- WetNatRich ~ bols(d90, intercept=FALSE)+bbs(d90, center=TRUE, df=1)+
bols(d365, intercept=FALSE)+bbs(d365, center=TRUE, df=1)+
bols(TSLW, intercept=FALSE)+bbs(TSLW, center=TRUE, df=1)+
bols(d1yrs_wet, intercept=FALSE)+bbs(d1yrs_wet, center=TRUE, df=1)+
bols(d3yrs_dry, intercept=FALSE)+bbs(d3yrs_dry, center=TRUE, df=1)+
bols(Inundated, intercept=FALSE)


# Run model on training set
Full.model <-gamboost(formulaB, family = MBztnegbin(), data = subset(datenSmall, WetNatRich > 0),control=boost_control(mstop=1000,trace=TRUE)) # originally 10000 but reduced down the mstop value for trial runs as it takes a while
# This churns out a couple of warnings that I don't fully understand regarding the linear effects - covariates should be (mean-) centered if inter


par(mfrow = c(4,3))    ## 3 plots in one device
plot(Full.model) 


g1 <- gamboost(satellites ~ bbs(width) + bbs(color),
  data = subset(CrabSatellites, satellites > 0), family = MBztnegbin())
set.seed(1)

# Carry out 5 fold cross validation to determine optimal stopping iteration - this seems to still be 10000 - increase cross validation for proper runs?
cv5f <- cv(model.weights(Full.model), type='subsampling', B=25)
cvm <- cvrisk(Full.model, folds=cv5f)
#plot(cvm)
st<-(mstop(cvm))
Full.model[st]
coef(Full.model)

predicted<-list() # creates a whole bunch of empty lists
predicted.insample <-list()
nuvec<-list()
null.model.coef<-list()
null.model.nu<-list()

# Predictions for out-of-bootstrap data
predictions_ofb<-predict(Full.model,newdata=testdata, type='response')
plot((predictions_ofb),(testdata$Wet_NativesPA))
abline(0,1)


######################################################################
# Plots

png(paste(image.dir,'H_index_fitted vs observed.png',sep=''), width=2000, height=1000, units="px", res=300)
par(mfrow=c(1,2))
par(mar=c(5,4,1,1))
plot(fitted(Full.model),datenSmall$Wet_NativesPA, xlab='Fitted values', ylab='Observed values')
abline(0,1)
plot(predict(Full.model,newdata=testdata),testdata$Wet_Natives, xlab='Predicted values', ylab='Observed values')
abline(0,1)
dev.off ()



png(paste(image.dir,'H_index_marginal_plots.png',sep=''), width=2000, height=2000, units="px", res=300)
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
# Marginal functional estimates of boosted additive models for flood frequency, time since last flood and rainfall in 90 days prior to sampling
# rem that data has been centred so when looking at the plots its helpful to 'uncentre' the data

# plot using d30(accumulated rainfall in 30 days prior to sampling)
md30<-mean(log10(data.matrix.env.data$d30+1))
xmatLin <- extract(Full.model, which=5)
xmatSmooth <- extract(Full.model, which=6)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d30, df = 1, center = TRUE)` + xmatLin[[1]]*coef(Full.model)$`bols(d30, intercept = FALSE)`
plot(sort(datenSmall$d30+md30),yvalues[order(datenSmall$d30+md30)], type="l",xlab='log(d30+1)', ylab='f(log(d30+1))')
rug(sort(datenSmall$d30+md30))


# plot using d365(accumulated rainfall in 90 days prior to sampling)
md365<-mean(log10(data.matrix.env.data$d365+1))
xmatSmooth <- extract(Full.model, which=8)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(d365, df = 1, center = TRUE)` 
plot(sort(datenSmall$d365+md365),yvalues[order(datenSmall$d365+md365)], type="l",xlab='log(d365+1)', ylab='f(log(d365+1))')
rug(sort(datenSmall$d365+md365))


mFf<-mean(envdata_backup$Flood_frequency)
xmatLin <- extract(Full.model, which=9)
xmatSmooth <- extract(Full.model, which=10)
# the below line had to adapted from the paper to correct for the type of speech marks, the order or the components and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(Flood_frequency, df = 1, center = TRUE)` + xmatLin[[1]]*coef(Full.model)$`bols(Flood_frequency, intercept = FALSE)`
plot(sort(datenSmall$Flood_frequency+mFf),yvalues[order(datenSmall$Flood_frequency+mFf)], type="l",xlab='Flood frequency', ylab='f(Flood frequency)')
rug(sort(datenSmall$Flood_frequency+mFf))

#  plot using Time Since Last Flood (TSLW)
mTSLW<-mean(log10(envdata_backup$TSLW+1))
xmatSmooth <- extract(Full.model, which=12)
# the below line had to adapted from the paper to correct for the type of speech marks, the order and the spacing around the equals signs: see coef(Full.model)
yvalues=xmatSmooth[[1]]%*%coef(Full.model)$`bbs(TSLW, df = 1, center = TRUE)`
plot(sort(datenSmall$TSLW+mTSLW),yvalues[order(datenSmall$TSLW+mTSLW)], type="l",xlab='log(TSLW+1)', ylab='f(log(TSLW+1))')
rug(sort(datenSmall$TSLW+mTSLW))



dev.off()


par(mfrow = c(3,4))
plot(Final.model, which = "bbs")
# plot just linear effects
par(mfrow = c(3,3))
plot(Full.model, which = "bols")
# Plot spatial effects



Calculate pseudo R2 using Nagelkerke

library(lme4)

m1 <-lmer(H_index~1+(1|Site.ID),data=datenSmall) # null model
#m1 <-lm(H_index~1,data=datenSmall) # null model

logLik.lmer.intercept=logLik(m1)
logLik.lmer.full=logLik(Full.model)
N.lmer.full=nrow(datenSmall)
coxsnell<- 1 - exp((logLik.lmer.intercept - logLik.lmer.full) * (2/N.lmer.full))
 nagelkerke<- coxsnell / (1 - exp(logLik.lmer.intercept * (2/N.lmer.full)))
 


R2gauss<- function(y,model){
    moy<-mean(y)
    N<- length(y)
    p<-length(coef(model))-1
    SSres<- sum((y-predict(model))^2)
    SStot<-sum((y-moy)^2)
    R2<-1-(SSres/SStot)
    Rajust<-1-(((1-R2)*(N-1))/(N-p-1))
    return(data.frame(R2,Rajust,SSres,SStot))
}