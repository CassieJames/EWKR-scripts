##############################################################################
#### Model richness and abundance of wetland and dryland species 
library(corrplot)
library(randomForest)
library(vegan)
library(forestFloor)


data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.env.data=read.csv("Hattah wetlands response by metrics.csv") # load data
image.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Plots/"

envdata=data.matrix.env.data
envdata$TSLW[envdata$Inundated==TRUE] <- 1 # If site is currently inundated overwrite TSLW and make it 1 day
envdata=envdata[!is.na(envdata$TSLW),] # remove sites where TSLW is beyond record
envdata$d3Mon_meandepth[is.na(envdata$d3Mon_meandepth)] <-0 # is the site has not flooded in a year make the mean depth 0
envdata$d1yrs_meandepth[is.na(envdata$d1yrs_meandepth)] <-0 # is the site has not flooded in a year make the mean depth 0
envdata$d3yrs_meandepth[is.na(envdata$d3yrs_meandepth)] <-0 # is the site has not flooded in a year make the mean depth 0


envdata.short <- envdata[,c("d365","d90", "TSLW", "MRI_length","MRI_meandepth","d1yrs_wet","d1yrs_meandepth","d3yrs_wet","d3yrs_meandepth", "Freq_d1","Freq_d3","MeanTemp90","MaxTemp90","MinTemp90")]


envdata.short.matrix<-as.matrix(as.data.frame(envdata.short)) 
M <- cor(envdata.short.matrix)
corrplot(M,type = "upper", tl.col = "black", tl.srt = 45) 

daten=envdata


Wet_Natives.mtry<-tuneRF(envdata.short, envdata$Wet_Natives, ntreeTry=1000, plot=TRUE)

Wet_Natives.rf<-randomForest(envdata$Wet_Natives ~ ., data=envdata.short, localImp=TRUE, 
      proximity=TRUE, importance=TRUE, ntree=1000, mtry=4, oob.prox=TRUE, 
      keep.forest=TRUE)

varImpPlot(Wet_Natives.rf, sort=TRUE, type=1)


par(mfrow=c(3,3))
par(mar=c(4,4,1,1))
partialPlot(Wet_Natives.rf,envdata.short, MRI_TSLW, main=" ", ylab="Native wetland species abundance")
partialPlot(Wet_Natives.rf,envdata.short, MeanTemp90, main=" ", ylab="Native wetland species abundance")
partialPlot(Wet_Natives.rf,envdata.short, MRI_meandepth, main=" ", ylab="Native wetland species abundance")