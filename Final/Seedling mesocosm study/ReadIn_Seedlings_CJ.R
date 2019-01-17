library(dplyr); library(tidyr)

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Mesocosm study/"; setwd(data.dir)

Seedlings <- read.csv('All_Mesocosm_Data.csv')

Seedlings$Mortality <- vector(mode = 'numeric', length = dim(Seedlings)[1])
Seedlings$RootsAbove <- vector(mode = 'numeric', length = dim(Seedlings)[1])
Seedlings$Mortality[Seedlings$`mortality` == 'l'] <- 0; Seedlings$Mortality[Seedlings$`mortality` == 'd'] <- 1
Seedlings$RootsAbove[Seedlings$`rootsAbove` == 'y'] <- 1; Seedlings$RootsAbove[Seedlings$`rootsAbove` == 'n'] <- 0
Seedlings$mortality <- NULL
Seedlings$rootsAbove <- NULL

# Subset only phase 2 of data for time being
Seedlings <- Seedlings[Seedlings$Week ==23, ]

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

# Initialise frames of summary stats by group:
SeedlingsLive %>%
  group_by(Treat, Species) %>%
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
            MassABSE = sd(MassAB, na.rm = TRUE)/sqrt(n()) ) -> SLstats

# Save data frames in named list:
SeedFrames <- list('Seedlings' = Seedlings, 'SeedlingsLive' = SeedlingsLive)
save(SeedFrames, file = 'Data_SeedlingsP2.RData')
