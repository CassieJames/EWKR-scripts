# Script to undertake basic MDS
# Written by  C.S.James 
# GNU General Public License .. feel free to use / distribute ... no warranties
# 29th July 2016
library(vegan)
library(labdsv)
library(RColorBrewer)

data.dir="C:/Users/jc246980/Documents/Documents (2)/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix=read.csv("Spp_Env_HTH_FP_Dec 2017.csv") # load data
Output=data.matrix[,c(3:267)] # need to check these numbers and change code eventually - just a cheat method for the time being
rownames(Output)=data.matrix$Row.names
remove=c("TR4SOT3_08", "TR5OFT3_14","TR5OFT4_14") # extreme outliers
Output=Output[!rownames(Output) %in% remove, ]
remove=c("Inundated.x", "Leaf litter >50%") # columns that are not live vegetation
Output=Output[,!colnames(Output) %in% remove ]
Output<-dropspc(Output, 2)   # remove species that occur in = < two sites/date combos
Output=Output[rowSums(Output!= 0) > 0,]	# remove sites with no records 

result<-metaMDS((Output), distance="bray", autotransform=F, k=2)
site.scores=scores(result, display="sites") 
tdata=as.data.frame(site.scores)
tdata$groups<- rownames(tdata)
tdata$groups[grepl("_08",tdata$groups)] <- "Y2008"
tdata$groups[grepl("_09",tdata$groups)] <- "Y2009"
tdata$groups[grepl("_10",tdata$groups)] <- "Y2010"
tdata$groups[grepl("_11",tdata$groups)] <- "Y2011"
tdata$groups[grepl("_12",tdata$groups)] <- "Y2012"
tdata$groups[grepl("_13",tdata$groups)] <- "Y2013"
tdata$groups[grepl("_14",tdata$groups)] <- "Y2014"
tdata$groups[grepl("_16",tdata$groups)] <- "Y2016"

tdata$groups=as.factor(tdata$groups)
colvec=colorRampPalette(c("lightblue","black"))(8)
plot(result, type = "n")
points(result,pch= 19, col=colvec[tdata$groups], cex=1.2)


groupz <- unique(tdata$groups)
for(i in seq(groupz)) {
ordisegments(result, tdata$groups, col=colvec[i], display = "sites", show.groups=groupz[i])
} 


# Mds plot with flood frequency plotted

site.scores=scores(result, display="sites") 
tdata=as.data.frame(site.scores)
tdata$groups<- rownames(tdata)
tdata$groups[grepl("OF",tdata$groups)] <- "OFTEN"
tdata$groups[grepl("SO",tdata$groups)] <- "SOMETIMES"
tdata$groups[grepl("RA",tdata$groups)] <- "RARELY"


tdata$groups=as.factor(tdata$groups)
colvec=c("Blue", "purple", "red")
plot(result, type = "n")
points(result,pch= 19, col=colvec[tdata$groups], cex=1.2)


groupz <- unique(tdata$groups)
for(i in seq(groupz)) {
ordisegments(result, tdata$groups, col=colvec[i], display = "sites", show.groups=groupz[i])
} 















