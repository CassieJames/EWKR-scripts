# script to plot basic mds and display by date and 
library(vegan)
library(labdsv)
library(RColorBrewer)

Output.com<-dropspc(Output, 2)   # remove species that occur in = > two sites/date combos
Output.com=Output.com[rowSums(Output.com != 0) > 0,]	# remove sites with no records

remove=c("TR4SOT3_08") # extreme outlier
Output.com=Output.com[!rownames(Output.com) %in% remove, ]

result<-metaMDS(log(Output.com+1), distance="bray", autotransform=F, k=2)
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


# With environmental data

Output.com<-dropspc(Output, 2)   # remove species that occur in = > two sites/date combos

ef1 <- envfit(result, envdata, permu = 999, na.rm=TRUE) 