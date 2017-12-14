library(dplyr)
library(gridExtra)

data.dir = "C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/SeedlingMesocosm/"; setwd (data.dir)
mydata=data.frame(read.csv("Data23Weeks.csv"))
mydata=mydata[which(mydata$Week==23),]

mydata$Treatment <-as.factor(mydata$Treatment)
mydata$Tank <-as.factor(mydata$Tank)
mydata$Species <-as.factor(mydata$Species)
mydata$Ratio <-as.numeric(mydata$Ratio)



Mydata_summary <- mydata %>% # the names of the new data frame and the data frame to be summarised
  group_by(Species,Treatment) %>%   # the grouping variable
  summarise(mean_SeedHt = mean(SeedlingHt),  # calculates the mean of each group
            sd_SeedHt = sd(SeedlingHt), # calculates the standard deviation of each group
            n_SeedHt = n(),  # calculates the sample size per group
            SE_SeedHt = sd(SeedlingHt)/sqrt(n()),
			mean_LeavesNo = mean(LeavesNo),  # calculates the mean of each group
            sd_LeavesNo = sd(LeavesNo), # calculates the standard deviation of each group
            n_LeavesNo = n(),  # calculates the sample size per group
            SE_LeavesNo = sd(LeavesNo)/sqrt(n()),
			mean_RootDepth = mean(RootDepth),  # calculates the mean of each group
            sd_RootDepth = sd(RootDepth), # calculates the standard deviation of each group
            n_RootDepth = n(),  # calculates the sample size per group
            SE_RootDepth = sd(RootDepth)/sqrt(n()),
			mean_Coppicing = mean(Coppicing),  # calculates the mean of each group
            sd_Coppicing = sd(Coppicing), # calculates the standard deviation of each group
            n_Coppicing = n(),  # calculates the sample size per group
            SE_Coppicing = sd(Coppicing)/sqrt(n()),
			Count_Dead = sum(ifelse(Mortality=='D',1,0)),
			Count_Alive = sum(ifelse(Mortality=='L',1,0)),
			Count_AirRoots = sum(ifelse(RootsAbove=='Y',1,0)),
			Count_NoAirRoots = sum(ifelse(RootsAbove=='N',1,0)),
			mean_AboveG = mean(AboveG),  # calculates the mean of each group
            sd_AboveG = sd(AboveG), # calculates the standard deviation of each group
            n_AboveG = n(),  # calculates the sample size per group
            SE_AboveG = sd(AboveG)/sqrt(n()),
			mean_BelowG = mean(BelowG),  # calculates the mean of each group
            sd_BelowG = sd(BelowG), # calculates the standard deviation of each group
            n_BelowG = n(),  # calculates the sample size per group
            SE_BelowG = sd(BelowG)/sqrt(n()),
			mean_ratio = mean(AboveG/BelowG,  na.rm=TRUE),  # calculates the mean of each group
            sd_ratio = sd(AboveG/BelowG,  na.rm=TRUE), # calculates the standard deviation of each group
            n_ratio = n(),  # calculates the sample size per group
            SE_ratio = sd(AboveG/BelowG,  na.rm=TRUE)/sqrt(n())
			
			) # calculates the standard error of each group

plt1<-ggplot(Mydata_summary, aes(x = Treatment, y = mean_SeedHt))+ ylab("Seedling height(cm)")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)
limits <- aes(ymax = mean_SeedHt +  SE_SeedHt,
              ymin = mean_SeedHt - SE_SeedHt)
plt1<-plt1 + geom_errorbar(limits, position = position_dodge(0.9),width = 0.1)

plt2<-ggplot(Mydata_summary, aes(x = Treatment, y = mean_LeavesNo))+ ylab("Leaf No.")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)
limits <- aes(ymax = mean_LeavesNo +  SE_LeavesNo,
              ymin = mean_LeavesNo - SE_LeavesNo)
plt2<-plt2 + geom_errorbar(limits, position = position_dodge(0.9),width = 0.1)

plt3<-ggplot(Mydata_summary, aes(x = Treatment, y = Count_Dead))+ ylab("Dead No.")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)


plt4<-ggplot(Mydata_summary, aes(x = Treatment, y = mean_RootDepth))+ ylab("Root depth(cm)")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)
limits <- aes(ymax = mean_RootDepth +  SE_RootDepth,
              ymin = mean_RootDepth - SE_RootDepth)
plt4<-plt4 + geom_errorbar(limits, position = position_dodge(0.9),width = 0.1)


plt5<-ggplot(Mydata_summary, aes(x = Treatment, y = mean_Coppicing))+ ylab("Coppicing No.")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)
limits <- aes(ymax = mean_Coppicing +  SE_Coppicing,
              ymin = mean_Coppicing - SE_Coppicing)
plt5<-plt5 + geom_errorbar(limits, position = position_dodge(0.9),width = 0.1)


plt6<-ggplot(Mydata_summary, aes(x = Treatment, y = Count_AirRoots))+ ylab("Aerial Root No.")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)


plt7<-ggplot(Mydata_summary, aes(x = Treatment, y = mean_AboveG))+ ylab("Above ground biomass (g)")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)
limits <- aes(ymax = mean_AboveG +  SE_AboveG,
              ymin = mean_AboveG - SE_AboveG)
plt7<-plt7 + geom_errorbar(limits, position = position_dodge(0.9),width = 0.1)



plt8<-ggplot(Mydata_summary, aes(x = Treatment, y = mean_BelowG))+ ylab("Below ground biomass (g)")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)
limits <- aes(ymax = mean_BelowG +  SE_BelowG,
              ymin = mean_BelowG - SE_BelowG)
plt8<-plt8 + geom_errorbar(limits, position = position_dodge(0.9),width = 0.1)



plt9<-ggplot(Mydata_summary, aes(x = Treatment, y = mean_ratio))+ ylab("Above/below ground ratio")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)
limits <- aes(ymax = mean_ratio +  SE_ratio,
              ymin = mean_ratio - SE_ratio)
plt9<-plt9 + geom_errorbar(limits, position = position_dodge(0.9),width = 0.1)



png(paste(data.dir,"measurements23Weeks.png",sep=''),width=12, height=30, units='cm', res=300, pointsize=15, bg='white')
grid.arrange(plt1,plt2,plt4,plt5,ncol=1)
dev.off()

png(paste(data.dir,"measurements23Weeks_mortality and Coppicing.png",sep=''),width=12, height=15, units='cm', res=300, pointsize=15, bg='white')
grid.arrange(plt3,plt6,ncol=1)
dev.off()

png(paste(data.dir,"measurements23Weeks_Biomass.png",sep=''),width=12, height=22.5, units='cm', res=300, pointsize=15, bg='white')
grid.arrange(plt7,plt8,plt9,ncol=1)
dev.off()

####################################################################################################################################

### Box plots

plt1<-ggplot(mydata, aes(x = Treatment, y = SeedlingHt))+ ylab("Seedling height(cm)")+geom_boxplot(fill='grey')+facet_wrap(~Species,ncol=3)
plt2<-ggplot(mydata, aes(x = Treatment, y = LeavesNo))+ ylab("Leaf No.")+geom_boxplot(fill='grey')+facet_wrap(~Species,ncol=3)
plt4<-ggplot(mydata, aes(x = Treatment, y = RootDepth))+ ylab("Root depth(cm)")+geom_boxplot(fill='grey')+facet_wrap(~Species,ncol=3)
plt5<-ggplot(mydata, aes(x = Treatment, y = Coppicing))+ ylab("Coppicing No.")+geom_boxplot(fill='grey')+facet_wrap(~Species,ncol=3)
plt7<-ggplot(mydata, aes(x = Treatment, y = AboveG))+ ylab("Above ground biomass (g)")+geom_boxplot(fill='grey')+facet_wrap(~Species,ncol=3)
plt8<-ggplot(mydata, aes(x = Treatment, y = BelowG))+ ylab("Below ground biomass (g)")+geom_boxplot(fill='grey')+facet_wrap(~Species,ncol=3)
plt9<-ggplot(mydata, aes(x = Treatment, y = Ratio))+ ylab("Above/below ground ratio")+geom_boxplot(fill='grey')+facet_wrap(~Species,ncol=3)+ coord_cartesian(ylim = c(0, 10))

plt3<-ggplot(Mydata_summary, aes(x = Treatment, y = Count_Dead))+ ylab("Dead No.")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)
plt6<-ggplot(Mydata_summary, aes(x = Treatment, y = Count_AirRoots))+ ylab("Aerial Root No.")+geom_bar(fill='grey',stat = "identity")+facet_wrap(~Species,ncol=3)

png(paste(data.dir,"23Weeks_measurements_boxplots.png",sep=''),width=15, height=30, units='cm', res=300, pointsize=15, bg='white')
grid.arrange(plt1,plt2,plt4,plt5,ncol=1)
dev.off()

png(paste(data.dir,"23Weeks_mortality and Coppicing.png",sep=''),width=15, height=15, units='cm', res=300, pointsize=15, bg='white')
grid.arrange(plt3,plt6,ncol=1)
dev.off()

png(paste(data.dir,"23Weeks_Biomass_boxplots.png",sep=''),width=15, height=22.5, units='cm', res=300, pointsize=15, bg='white')
grid.arrange(plt7,plt8,plt9,ncol=1)
dev.off()



####################################################################################################################################
#### Analysis of survival

library(lattice)
library(lme4)
library(alr3)
library(nlme)

survival.glm <- glmer(Mortality ~ Treatment * Species+(1|Tank), family=binomial,data=mydata)

survival.glm <- glmer(Mortality ~ Tank + Species +(1|Tank), family=binomial,data=mydata)


library(blmeco) 
dispersion_glmer(survival.glm)

par(mfrow=c(2,2))
qqnorm(resid(survival.glm), main="normal qq-plot, residuals")
qq
line(resid(survival.glm))

mydata_alive=mydata[which(mydata$Mortality==1),]


SeedHt.glm <- lmer(SeedlingHt ~ Treatment * Species+(1|Tank),data=mydata)