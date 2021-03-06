library(MASS)
library(lme4)
library(nlme)
library(lmtest)
library(emmeans)

data.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Mesocosm study/"; setwd(data.dir)
load('Data_SeedlingsP2.RData')
SL <- SeedFrames$SeedlingsLive 


nLeaves.lme1 <- lme(nLeaves ~ Species * Treat, data = SL, random = ~ 1 | Tank)
#nLeaves.glme1 <- glmmPQL(nLeaves ~ Species * Treat, data = SL, random = ~ 1 | Tank, family = poisson)
nLeaves.glmeNull<- glmer(nLeaves ~ 1 | Tank, data = SL, family  = poisson)
nLeaves.glmeSpecies <- glmer(nLeaves ~ Species + (1 | Tank), data = SL, family  = poisson)
nLeaves.glmeTreat <- glmer(nLeaves ~ Treat + (1 | Tank), data = SL, family  = poisson)
nLeaves.glmeNoInt <- glmer(nLeaves ~ Species + Treat + (1 | Tank), data = SL, family  = poisson)
nLeaves.glmeFull<- glmer(nLeaves ~ Species * Treat + (1 | Tank), data = SL, family  = poisson)
nLeavesAICs <- AIC(nLeaves.lme1, nLeaves.glmeNull, nLeaves.glmeTreat, nLeaves.glmeSpecies, nLeaves.glmeNoInt, nLeaves.glmeFull)

anova(nLeaves.glmeFull, nLeaves.glmeNoInt, test = 'Chi') # interaction is significant

lArea.lme1 <- lme(LeafArea ~ Species * Treat, data = SL, random = ~ 1 | Tank)
lArea.lme2  <- lme(LeafArea ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
lArea.lme2.1  <- lme(LeafArea ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
lArea.lme3  <- lme(LeafArea ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
lArea.lme4  <- lme(LeafArea ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
lArea.lme5  <- lme(LeafArea ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
lArea.lme6  <- lme(LeafArea ~ Species, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
lArea.lme7  <- lme(LeafArea ~ Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
lAreaAICs <- AIC(lArea.lme1, lArea.lme2, lArea.lme3, lArea.lme4, lArea.lme5, lArea.lme6, lArea.lme7)

lArea.lme2  <- lme(LeafArea ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat), method='ML')
lArea.lme5  <- lme(LeafArea ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat),method='ML')
anova(lArea.lme2, lArea.lme5) #L=98.13 # leaf area but no interaction

AverageLA.lme1 <- lme(AverageLA ~ Species * Treat, data = SL, random = ~ 1 | Tank)
AverageLA.lme2  <- lme(AverageLA ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
AverageLA.lme2.1  <- lme(AverageLA ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
AverageLA.lme3  <- lme(AverageLA ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
AverageLA.lme3.1  <- lme(AverageLA ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
AverageLA.lme4  <- lme(AverageLA ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
AverageLA.lme5  <- lme(AverageLA ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
AverageLA.lme6  <- lme(AverageLA ~ Species, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
AverageLA.lme7  <- lme(AverageLA ~ Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
AverageLAAICs <- AIC(AverageLA.lme1, AverageLA.lme2, AverageLA.lme3,AverageLA.lme3.1, AverageLA.lme4, AverageLA.lme5, AverageLA.lme6, AverageLA.lme7)

AverageLA.lme3  <- lme(AverageLA ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species), method='ML')
AverageLA.lme3.1  <- lme(AverageLA ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species), method='ML')
anova(AverageLA.lme3, AverageLA.lme3.1) 


emmeans(lArea.lme5, list(pairwise ~ Species), adjust = "tukey") # RRG differed from other two species
emmeans(lArea.lme5, list(pairwise ~ Treat), adjust = "tukey") # 1 is diff from 3,4 , 2 is diff from 3 and 5 is diff from 1,2,3,4

Height.lme1 <- lme(Height ~ Species * Treat, data = SL, random = ~ 1 | Tank)
Height.lme2  <- lme(Height ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
Height.lme3  <- lme(Height ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
Height.lme4  <- lme(Height ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
Height.lme5  <- lme(Height ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
Height.lme6  <- lme(Height ~ Species, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
Height.lme7  <- lme(Height ~ Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
HeightAICs <- AIC(Height.lme1, Height.lme2, Height.lme3, Height.lme4, Height.lme5, Height.lme6, Height.lme7)

Height.lme4  <- lme(Height ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat),method='ML')
Height.lme5  <- lme(Height ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat),method='ML')
anova(Height.lme4, Height.lme5) # L=56.73 and interaction not significant
emmeans(Height.lme5, list(pairwise ~ Species), adjust = "tukey") #species all different
emmeans(Height.lme5, list(pairwise ~ Treat), adjust = "tukey") # only NOT sig differences are 1 not diff from 2 and, 3 not diff from 4

RootLength.lme1 <- lme(RootLength ~ Species * Treat, data = SL, random = ~ 1 | Tank)
RootLength.lme2  <- lme(RootLength ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
RootLength.lme3  <- lme(RootLength ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
RootLength.lme4  <- lme(RootLength ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
RootLength.lme5  <- lme(RootLength ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
RootLength.lme6  <- lme(RootLength ~ Species, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
RootLength.lme7  <- lme(RootLength ~ Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
RootLengthAICs <- AIC(RootLength.lme1, RootLength.lme2, RootLength.lme3, RootLength.lme4, RootLength.lme5, RootLength.lme6, RootLength.lme7)

RootLength.lme4  <- lme(RootLength ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat),method='ML')
RootLength.lme5  <- lme(RootLength ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat),method='ML')
anova(RootLength.lme4, RootLength.lme5) # L=58.09 and interaction not significant
emmeans(RootLength.lme5, list(pairwise ~ Species), adjust = "tukey") ## RRG greater root lenth than BB and C
emmeans(RootLength.lme5, list(pairwise ~ Treat), adjust = "tukey") # 2 NOT diff from 3 and 4 and, 3 and 4 NOT diff

Coppice.lme1 <- lme(Coppice ~ Species * Treat, data = SL, random = ~ 1 | Tank)
#Coppice.glme1 <- glmmPQL(Coppice ~ Species * Treat, data = SL, random = ~ 1 | Tank, family = poisson)
Coppice.glmeNull<- glmer(Coppice ~ 1 | Tank, data = SL, family  = poisson)
Coppice.glmeSpecies <- glmer(Coppice ~ Species + (1 | Tank), data = SL, family  = poisson)
Coppice.glmeTreat <- glmer(Coppice ~ Treat + (1 | Tank), data = SL, family  = poisson)
Coppice.glmeNoInt <- glmer(Coppice ~ Species + Treat + (1 | Tank), data = SL, family  = poisson)
Coppice.glmeFull<- glmer(Coppice ~ Species * Treat + (1 | Tank), data = SL, family  = poisson)
CoppiceAICs <- AIC(Coppice.lme1, Coppice.glmeNull, Coppice.glmeTreat, Coppice.glmeSpecies, Coppice.glmeNoInt, Coppice.glmeFull)

anova(Coppice.glmeFull, Coppice.glmeNoInt, test = 'chi') # Interaction significant

MassAbove.lme1 <- lme(MassAbove ~ Species * Treat, data = SL, random = ~ 1 | Tank)
MassAbove.lme2  <- lme(MassAbove ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
MassAbove.lme3  <- lme(MassAbove ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
MassAbove.lme4  <- lme(MassAbove ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
MassAbove.lme5  <- lme(MassAbove ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat)) # Convergence probs with full-factor varIdent, so opt for next-down varIdent explaining most variance
MassAbove.lme6  <- lme(MassAbove ~ Species, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
MassAbove.lme7  <- lme(MassAbove ~ Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
MassAboveAICs <- AIC(MassAbove.lme1, MassAbove.lme2, MassAbove.lme3, MassAbove.lme4, MassAbove.lme5, MassAbove.lme6, MassAbove.lme7)

MassAbove.lme4  <- lme(MassAbove ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat),method='ML')
MassAbove.lme5  <- lme(MassAbove ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat),method='ML') 

anova(MassAbove.lme4, MassAbove.lme5)#L = 61.15 and no significant interaction
emmeans(MassAbove.lme5, list(pairwise ~ Species), adjust = "tukey") # species all different
emmeans(MassAbove.lme5, list(pairwise ~ Treat), adjust = "tukey") # 1 and 2 NOT diff and, 3 and 4 NOT diff

MassBelow.lme1 <- lme(MassBelow ~ Species * Treat, data = SL, random = ~ 1 | Tank)
MassBelow.lme2  <- lme(MassBelow ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
MassBelow.lme3  <- lme(MassBelow ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
MassBelow.lme4  <- lme(MassBelow ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
MassBelow.lme5  <- lme(MassBelow ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat))
MassBelow.lme6  <- lme(MassBelow ~ Species, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
MassBelow.lme7  <- lme(MassBelow ~ Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
MassBelowAICs <- AIC(MassBelow.lme1, MassBelow.lme2, MassBelow.lme3, MassBelow.lme4, MassBelow.lme5, MassBelow.lme6, MassBelow.lme7)

MassBelow.lme4  <- lme(MassBelow ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat),method='ML')
MassBelow.lme5  <- lme(MassBelow ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species*Treat),method='ML')
anova(MassBelow.lme4, MassBelow.lme5) #L=49.84 and no significant interaction
emmeans(MassBelow.lme5, list(pairwise ~ Species), adjust = "tukey") # RGR different
emmeans(MassBelow.lme5, list(pairwise ~ Treat), adjust = "tukey") # 1 and 2 NOT diff and, 3 and 4 NOT diff

MassAB.lme1 <- lme(MassAB ~ Species * Treat, data = SL, random = ~ 1 | Tank)
MassAB.lme2  <- lme(MassAB ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
MassAB.lme3  <- lme(MassAB ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
MassAB.lme4  <- lme(MassAB ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat)) # Convergence probs so simplification
MassAB.lme5  <- lme(MassAB ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
MassAB.lme6  <- lme(MassAB ~ Species, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Species))
MassAB.lme7  <- lme(MassAB ~ Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat))
MassABAICs <- AIC(MassAB.lme1, MassAB.lme2, MassAB.lme3, MassAB.lme4, MassAB.lme5, MassAB.lme6, MassAB.lme7)

MassAB.lme4  <- lme(MassAB ~ Species * Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat),method='ML')
MassAB.lme5  <- lme(MassAB ~ Species + Treat, data = SL, random = ~ 1 | Tank, weights = varIdent(form =~ 1 | Treat),method='ML')

anova(MassAB.lme4, MassAB.lme5) # no significant interaction

emmeans(MassAB.lme5, list(pairwise ~ Species), adjust = "tukey")
 # significant differences between BB and other two species, other two species NS different
 
 emmeans(MassAB.lme5, list(pairwise ~ Treat), adjust = "tukey")
 # significant differences between treatment 4 and 5 (DDDD and FFFF) 
