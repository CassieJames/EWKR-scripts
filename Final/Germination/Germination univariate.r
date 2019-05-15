#########################################################################################################################
#### Subset data into the balanced component. subset 1 is FF C2, C3 and C4 for Veg IS and NWW at locations LM, MQ and NL

library(lme4)
library(emmeans)
library(knitr)
library(car)
library(MASS)
library(vcd)
library(AER)
library(pscl)
library(rcompanion)
library(betareg)


data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Germination results/"; setwd (data.dir)
mytada=read.csv("Germination univariate summary.csv") 
subset1 <- mytada[!mytada$Location=="Middle Murray",]
subset1 <- subset1[!subset1$Flow_cat=="C1",]
subset1 <- subset1[!subset1$Veg=="IW",]

# set the contrasts for Anova

contrasts(subset1$Location) <- contr.sum 
contrasts(subset1$Flow_cat) <- contr.sum
contrasts(subset1$Veg) <- contr.sum

# in overdispersed models, the qq trend will deviate substantially from a straight line
# non-linear models will display trends in the residuals
# Goodness of fit 1 - pchisq(dat.glmL$deviance, dat.glmL$df.resid)
# 1 - pchisq(dat.resid, dat.glmL$df.resid) # Pearson's χ2 residuals - explores whether there are any significant patterns remaining in the residuals

#### Seedling abundance
# The paper of Ohara and Kotze (2010) recommend not logging abundance data in most circumstances however the data is highly overdispersed (so poisson not suitable)

r <- c(mean(subset1$Abund), var(subset1$Abund))
c(mean=r[1], var=r[2], ratio=r[2]/r[1]) # Note that variance is huge relative to the mean so poisson unlikely to be a good fit

# poisson
abund.glm.pois <- glm(Abund ~ Location * Flow_cat*Veg, data = subset1, family='poisson')
1-pchisq(abund.glm.pois$deviance,abund.glm.pois$df.residual) # p = zero which suggests terrible fit

#quasipoisson
abund.glmquaspois <- glm(Abund ~ Location * Flow_cat*Veg, data = subset1, family=quasipoisson(link=log))
res <- residuals(abund.glmquaspois, type="deviance")
plot(log(predict(abund.glmquaspois)), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)
1-pchisq(abund.glmquaspois$deviance,abund.glmquaspois$df.residual) # 

#negative binomial	
abund.glmnb <- glm.nb(Abund ~ Location * Flow_cat*Veg, data = subset1)
res <- residuals(abund.glmnb, type="deviance")
plot(log(predict(abund.glmnb)), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)	
1-pchisq(abund.glmnb$deviance,abund.glmnb$df.residual) # p = 0.028 is significant but a much better fit compared with poisson or quasipoisson

abund.glmnb.1 = update(abund.glmnb, . ~ . - Location:Flow_cat:Veg)
abund.glmnb.2 = update(abund.glmnb.1, . ~ . - Flow_cat:Veg)
abund.glmnb.3 = update(abund.glmnb.2, . ~ . - Location:Veg)

capture.output(anova(abund.glmnb.3,test="Chisq"), file="Species_abundance.txt") 

emms1 <- emmeans(abund.glmnb,pairwise ~ Location | Flow_cat, type="response")
emms2 <- emmeans(abund.glmnb,pairwise ~ Flow_cat | Location, type="response") # low category conditioned by location
emms3 <- emmeans(abund.glmnb,pairwise ~ Veg,type="response")

#### Species richness
par(mfrow = c(2, 2), mar = c(3, 3, 1, 1))

nSpecies.glm.null <- glm(Rich ~ 1, data = subset1, family  = poisson)
nSpecies.glm.full <- glm(Rich ~ Location * Flow_cat *Veg, data = subset1, family  = poisson)
#nSpecies.glm.full <- glm(Rich ~ Location * Flow_cat *Veg, data = subset1, family=quasipoisson(link=log))

1 - pchisq(nSpecies.glm.full$deviance, nSpecies.glm.full$df.resid)

dat.resid <- sum(resid(nSpecies.glm.full, type = "pearson")^2)
1 - pchisq(dat.resid, nSpecies.glm.full$df.resid)

dispersiontest(nSpecies.glm.full)
influencePlot(nSpecies.glm.full)
res <- residuals(nSpecies.glm.full, type="deviance")
plot(log(predict(nSpecies.glm.full)), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)

anova(nSpecies.glm.full,test="Chisq")
capture.output(anova(nSpecies.glm.full,test="Chisq"), file="Species_Richness.txt") 
emmeans(nSpecies.glm.full1, pairwise ~ Location,type="response")

#### Exotic proportion
# Data is very difficult to analyse as lots of zeros and then some high values and its not possible to parameterize a zero inflated model with the number of predictors.
# Kruskal Wallis test is for one way nonparametric tests only so not interactions and multiple factors

exotics = cbind(subset1$ExoticAbund, subset1$Abund-subset1$ExoticAbund)  
Exoticpa.glm.fullv2 <- glm(exotics ~ Location * Flow_cat*Veg, data = subset1, family=binomial(link="logit"))
Anova(Exoticpa.glm.fullv2,type="II",test="Wald")

dat.sim <- simulate(Exoticpa.glm.fullv2, n = 250)
par(mfrow = c(5, 4), mar = c(3, 3, 1, 1))
resid <- NULL
for (i in 1:nrow(dat.sim)) {
    e = ecdf(data.matrix(dat.sim[i, ] + runif(250, -0.5, 0.5)))
    plot(e, main = i, las = 1)
    resid[i] <- e(subset1$ExoticAbund/subset1$Abund[i] + runif(250, -0.5, 0.5))
}

plot(resid ~ fitted(Exoticpa.glm.fullv2))

dat.resid <- sum(resid(Exoticpa.glm.fullv2, type = "pearson")^2)
1 - pchisq(dat.resid, Exoticpa.glm.fullv2$df.resid)
1 - pchisq(Exoticpa.glm.fullv2$deviance, Exoticpa.glm.fullv2$df.resid)
# see http://www.flutterbys.com.au/stats/tut/tut10.5a.html
# Residuals appear to have very strong pattern and both the pearsons r2 and deviance tests are highly significant suggesting poor fit
# Think that the only option is to base exotics analysis on pa analysis

#### Native abundance
# The paper of Ohara and Kotze (2010) recommend not logging abundance data in most circumstances however the data is highly overdispersed and the negative binomial model 
NativeAbund.aov1 <- aov(log10(NativeAbund) ~ Location * Flow_cat*Veg, data = subset1)
op <-  par(mfrow = c(2, 2))
plot(NativeAbund.aov1)
par(op)
anova(NativeAbund.aov1) # no significant interactions just effect of flow cat and location
capture.output(anova(NativeAbund.aov1), file="Native_Abundance.txt") 
emmeans(NativeAbund.aov1, pairwise ~ Location)
emmeans(NativeAbund.aov1, pairwise~Flow_cat)
emmeans(NativeAbund.aov1, pairwise~Veg)

#### Native richness

nSpecies.glm.null <- glm(Rich ~ 1, data = subset1, family  = poisson)
nSpecies.glm.full <- glm(Rich ~ Location * Flow_cat *Veg, data = subset1, family  = poisson)

deviance(nSpecies.glm.full)/nSpecies.glm.full$df.residual
dispersiontest(nSpecies.glm.full)
influencePlot(nSpecies.glm.full)
res <- residuals(nSpecies.glm.full, type="deviance")
plot(log(predict(nSpecies.glm.full)), res)
abline(h=0, lty=2)
qqnorm(res)
qqline(res)

anova(nSpecies.glm.full,test="Chisq")
capture.output(anova(nSpecies.glm.full,test="Chisq"), file="Species_Richness.txt") 
emmeans(nSpecies.glm.full, pairwise ~ Location)


