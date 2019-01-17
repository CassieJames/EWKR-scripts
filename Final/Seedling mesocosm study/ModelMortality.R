library(nlme)
library(MASS)
library(lme4)

load('Data_Seedlings.RData')
Seedlings <- SeedFrames$Seedlings

mort.lme1 <- lme(Mortality ~ Species * Treat, data = Seedlings, random = ~ 1 | Tank)
#mort.glme1 <- glmmPQL(Mortality ~ Species * Treat, data = Seedlings, random = ~ 1 | Tank, family = binomial)
mort.glmeFull <- glmer(Mortality ~ Species * Treat + (1 | Tank), data = Seedlings, family = binomial)
mort.glmeNoInt <- glmer(Mortality ~ Species + Treat + (1 | Tank), data = Seedlings, family = binomial)
mort.glmeSpecies <- glmer(Mortality ~ Species + (1 | Tank), data = Seedlings, family = binomial)
mort.glmeTreat <- glmer(Mortality ~ Treat + (1 | Tank), data = Seedlings, family = binomial)
mort.glmeNull <- glmer(Mortality ~ 1 | Tank, data = Seedlings, family = binomial)

MortalityAICs <- AIC(mort.lme1, mort.glmeNull, mort.glmeSpecies, mort.glmeTreat, mort.glmeNoInt,  mort.glmeFull)

