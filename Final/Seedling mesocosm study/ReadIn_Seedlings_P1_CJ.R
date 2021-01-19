library(dplyr); library(tidyr)

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Mesocosm study/"; setwd(data.dir)

Seedlings <- read.csv('All_Mesocosm_Data.csv')


Calculate averages for week 3
Seedlingswk3 <- Seedlings[which(Seedlings$Week==3),]

Seedlingswk3$Mortality <- vector(mode = 'numeric', length = dim(Seedlingswk3)[1])
Seedlingswk3$RootsAbove <- vector(mode = 'numeric', length = dim(Seedlingswk3)[1])
Seedlingswk3$Mortality[Seedlingswk3$`mortality` == 'l'] <- 0; Seedlingswk3$Mortality[Seedlingswk3$`mortality` == 'd'] <- 1
SeedlingsLivewk3 <- Seedlingswk3[Seedlingswk3$Tank != 15 & Seedlingswk3$Mortality != 1, ]
SeedlingsLivewk3$MassAB <- SeedlingsLivewk3$MassAbove / SeedlingsLivewk3$MassBelow

# Initialise factors/nominals:
SeedlingsLivewk3$Species <- factor(SeedlingsLivewk3$Species)
SeedlingsLivewk3$Treat <- factor(SeedlingsLivewk3$Treat)
SeedlingsLivewk3$Treat_newP1 <- factor(SeedlingsLivewk3$Treat_newP1)
SeedlingsLivewk3$Tank <- factor(SeedlingsLivewk3$Tank)
SeedlingsLivewk3$Species <- factor(SeedlingsLivewk3$Species)
SeedlingsLivewk3$Treat <- factor(SeedlingsLivewk3$Treat)
SeedlingsLivewk3$Tank <- factor(SeedlingsLivewk3$Tank)
SeedlingsLivewk3 %>%
  group_by(Species) %>%
  summarise(nLeavesMean = mean(nLeaves, na.rm = TRUE),
            nLeavesSD = sd(nLeaves, na.rm = TRUE),
            nLeavesSE = sd(nLeaves, na.rm = TRUE)/sqrt(n()),
            CoppiceMean = mean(Coppice, na.rm = TRUE),
            CoppiceSD = sd(Coppice, na.rm = TRUE),
            CoppiceSE = sd(Coppice, na.rm = TRUE)/sqrt(n()),
            LeafAreaMean = mean(LeafArea, na.rm = TRUE),
            LeafAreaSD = sd(LeafArea, na.rm = TRUE),
            LeafAreaSE = sd(LeafArea, na.rm = TRUE)/sqrt(n()),
            HeightMean = mean(Height, na.rm = TRUE),
            HeightSD = sd(Height, na.rm = TRUE),
            HeightSE = sd(Height, na.rm = TRUE)/sqrt(n()),
            RootLengthMean = mean(RootLength, na.rm = TRUE),
            RootLengthSD = sd(RootLength, na.rm = TRUE),
            RootLengthSE = sd(RootLength, na.rm = TRUE)/sqrt(n()),
            MassAboveMean = mean(MassAbove, na.rm = TRUE),
            MassAboveSD = sd(MassAbove, na.rm = TRUE),
            MassAboveSE = sd(MassAbove, na.rm = TRUE)/sqrt(n()),
            MassBelowMean = mean(MassBelow, na.rm = TRUE),
            MassBelowSD = sd(MassAbove, na.rm = TRUE),
            MassBelowSE = sd(MassAbove, na.rm = TRUE)/sqrt(n()),
            MassABMean = mean(MassAB, na.rm = TRUE),
            MassABSD = sd(MassAB, na.rm = TRUE),
            MassABSE = sd(MassAB, na.rm = TRUE)/sqrt(n()), 
			MassTotalMean=mean(MassAbove+MassBelow,na.rm = TRUE)) -> SLstatswk3

SLstatswk3=as.data.frame(SLstatswk3)

Seedlings$Mortality <- vector(mode = 'numeric', length = dim(Seedlings)[1])
Seedlings$RootsAbove <- vector(mode = 'numeric', length = dim(Seedlings)[1])
Seedlings$Mortality[Seedlings$`mortality` == 'l'] <- 0; Seedlings$Mortality[Seedlings$`mortality` == 'd'] <- 1
Seedlings$RootsAbove[Seedlings$`rootsAbove` == 'y'] <- 1; Seedlings$RootsAbove[Seedlings$`rootsAbove` == 'n'] <- 0
Seedlings$mortality <- NULL
Seedlings$rootsAbove <- NULL

# Subset only phase 2 of data for time being
Seedlings <- Seedlings[Seedlings$Week ==13, ]

# Initialise another frame not containing rows of dead seedlings:
SeedlingsLive <- Seedlings[Seedlings$Tank != 15 & Seedlings$Mortality != 1, ]

# Add another response to SeedlingsLive:
SeedlingsLive$MassAB <- SeedlingsLive$MassAbove / SeedlingsLive$MassBelow

# Initialise factors/nominals:
Seedlings$Species <- factor(Seedlings$Species)
Seedlings$Treat <- factor(Seedlings$Treat)
Seedlings$Treat_newP1 <- factor(Seedlings$Treat_newP1)
Seedlings$Tank <- factor(Seedlings$Tank)
SeedlingsLive$Species <- factor(SeedlingsLive$Species)
SeedlingsLive$Treat <- factor(SeedlingsLive$Treat)
SeedlingsLive$Tank <- factor(SeedlingsLive$Tank)

# Determine growth rates using summary stats from week 3

metrics <- c('nLeaves','Coppice','LeafArea','Height','RootLength','MassAbove', 'MassBelow','MassAB', 'MassTotal')
tada = matrix(NA,nrow=nrow(SeedlingsLive),ncol=length(metrics)+5)#define the output matrix
colnames(tada)=c("Species","Treat","Treat_new","Treat_newP1","Tank",metrics)
tada=as.data.frame(tada)
tada$Species <- factor(SeedlingsLive$Species)
tada$Treat <- factor(SeedlingsLive$Treat)
tada$Tank <- factor(SeedlingsLive$Tank)
tada$Treat_new <- factor(SeedlingsLive$Treat_new)
tada$Treat_newP1 <- factor(SeedlingsLive$Treat_newP1)


Species <- factor(SeedlingsLive$Species)
metrics <- c('nLeaves','Coppice','LeafArea','Height','RootLength','MassAbove', 'MassBelow','MassAB', 'MassTotal')

for (s in Species) {
tdata <- SeedlingsLive[SeedlingsLive$Species ==s, ]
for (m in metrics) {
ttdata <- tdata[,colnames(tdata)==m]
wk3value=SLstatswk3[which(SLstatswk3$Species==s),colnames(SLstatswk3)==paste(m,"Mean",sep="")]
RGR=(ttdata-wk3value)
tada[grep(s,tada$Species),grep(m,colnames(tada))] <-RGR
}}

save(tada, file = 'Data_SeedlingsP1_RGR.RData')

# Initialise frames of summary stats by group:
tada %>%
  group_by(Treat_newP1, Species) %>%
  summarise(nLeavesMean = mean(nLeaves, na.rm = TRUE),
            nLeavesSD = sd(nLeaves, na.rm = TRUE),
            nLeavesSE = sd(nLeaves, na.rm = TRUE)/sqrt(n()),
            CoppiceMean = mean(Coppice, na.rm = TRUE),
            CoppiceSD = sd(Coppice, na.rm = TRUE),
            CoppiceSE = sd(Coppice, na.rm = TRUE)/sqrt(n()),
			LeafAreaMean = mean(LeafArea, na.rm = TRUE),
            LeafAreaSD = sd(LeafArea, na.rm = TRUE),
            LeafAreaSE = sd(LeafArea, na.rm = TRUE)/sqrt(n()),
            HeightMean = mean(Height, na.rm = TRUE),
            HeightSD = sd(Height, na.rm = TRUE),
            HeightSE = sd(Height, na.rm = TRUE)/sqrt(n()),
            RootLengthMean = mean(RootLength, na.rm = TRUE),
            RootLengthSD = sd(RootLength, na.rm = TRUE),
            RootLengthSE = sd(RootLength, na.rm = TRUE)/sqrt(n()),
            MassAboveMean = mean(MassAbove, na.rm = TRUE),
            MassAboveSD = sd(MassAbove, na.rm = TRUE),
            MassAboveSE = sd(MassAbove, na.rm = TRUE)/sqrt(n()),
            MassBelowMean = mean(MassBelow, na.rm = TRUE),
            MassBelowSD = sd(MassBelow, na.rm = TRUE),
            MassBelowSE = sd(MassBelow, na.rm = TRUE)/sqrt(n()),
			MassTotalMean = mean(MassTotal, na.rm = TRUE),
            MassTotalSD = sd(MassTotal, na.rm = TRUE),
            MassTotalSE = sd(MassTotal, na.rm = TRUE)/sqrt(n()),
            MassABMean = mean(MassAB, na.rm = TRUE),
            MassABSD = sd(MassAB, na.rm = TRUE),
            MassABSE = sd(MassAB, na.rm = TRUE)/sqrt(n()) ) -> SLstatsP1RGR




# Initialise frames of summary stats by group:
SeedlingsLive %>%
  group_by(Treat_newP1, Species) %>%
  summarise(nLeavesMean = mean(nLeaves, na.rm = TRUE),
            nLeavesSD = sd(nLeaves, na.rm = TRUE),
            nLeavesSE = sd(nLeaves, na.rm = TRUE)/sqrt(n()),
            CoppiceMean = mean(Coppice, na.rm = TRUE),
            CoppiceSD = sd(Coppice, na.rm = TRUE),
            CoppiceSE = sd(Coppice, na.rm = TRUE)/sqrt(n()),
            AverageLAMean= mean(AverageLA, na.rm = TRUE),
			AverageLASD= sd(AverageLA, na.rm = TRUE),
			AverageLASE= sd(AverageLA, na.rm = TRUE)/sqrt(n()),
			LeafAreaMean = mean(LeafArea, na.rm = TRUE),
            LeafAreaSD = sd(LeafArea, na.rm = TRUE),
            LeafAreaSE = sd(LeafArea, na.rm = TRUE)/sqrt(n()),
            HeightMean = mean(Height, na.rm = TRUE),
            HeightSD = sd(Height, na.rm = TRUE),
            HeightSE = sd(Height, na.rm = TRUE)/sqrt(n()),
            RootLengthMean = mean(RootLength, na.rm = TRUE),
            RootLengthSD = sd(RootLength, na.rm = TRUE),
            RootLengthSE = sd(RootLength, na.rm = TRUE)/sqrt(n()),
            MassAboveMean = mean(MassAbove, na.rm = TRUE),
            MassAboveSD = sd(MassAbove, na.rm = TRUE),
            MassAboveSE = sd(MassAbove, na.rm = TRUE)/sqrt(n()),
            MassBelowMean = mean(MassBelow, na.rm = TRUE),
            MassBelowSD = sd(MassAbove, na.rm = TRUE),
            MassBelowSE = sd(MassAbove, na.rm = TRUE)/sqrt(n()),
            MassABMean = mean(MassAB, na.rm = TRUE),
            MassABSD = sd(MassAB, na.rm = TRUE),
            MassABSE = sd(MassAB, na.rm = TRUE)/sqrt(n()) ) -> SLstats

# Save data frames in named list:
SeedFrames <- list('Seedlings' = Seedlings, 'SeedlingsLive' = SeedlingsLive)
save(SeedFrames, file = 'Data_SeedlingsP1.RData')
