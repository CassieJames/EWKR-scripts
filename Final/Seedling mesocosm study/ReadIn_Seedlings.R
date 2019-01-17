library(readxl); library(dplyr); library(tidyr)

Seedlings <- read_excel('Seedlings.xlsx', sheet = 'week 23_data all', 
                        col_types = c(rep('numeric',4), 'text', rep('numeric', 3), 'text', rep('numeric',2), 'text', rep('numeric',4)),
                        col_names = c('Pot','Week','Treat','Tank','Species','SoilMoisture','Height','nLeaves','mortality','RootLength',
                                      'Coppice','rootsAbove','MassAbove','MassBelow','MassTotal','LeafArea'), skip = 1)
Seedlings$Mortality <- vector(mode = 'numeric', length = dim(Seedlings)[1])
Seedlings$RootsAbove <- vector(mode = 'numeric', length = dim(Seedlings)[1])
Seedlings$Mortality[Seedlings$`mortality` == 'L'] <- 0; Seedlings$Mortality[Seedlings$`mortality` == 'D'] <- 1
Seedlings$RootsAbove[Seedlings$`rootsAbove` == 'Y'] <- 1; Seedlings$RootsAbove[Seedlings$`rootsAbove` == 'N'] <- 0
Seedlings$mortality <- NULL
Seedlings$rootsAbove <- NULL

# Initialise another frame not containing rows of dead seedlings:
SeedlingsLive <- Seedlings[Seedlings$Tank != 15 & Seedlings$Mortality != 1, ]

# Add another response to SeedlingsLive:
SeedlingsLive$MassAB <- SeedlingsLive$MassAbove / SeedlingsLive$MassBelow

# Initialise factors/nominals:
Seedlings$Species <- factor(Seedlings$Species)
Seedlings$Treat <- factor(Seedlings$Treat)
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
save(SeedFrames, file = 'Data_Seedlings.RData')
