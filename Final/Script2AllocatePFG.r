######################################################################################################################################
### Script used to allocated PFG and weed status to species list

library(vegan)
library(labdsv)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(BAT)

######################################################################################################################################
### Step 1 load data from different files (some is summarised to wetland by year and others to wetland with years merged already)

data.dir="C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd(data.dir)
image.dir="C:/Users/jc246980/Documents/MD Vegetation/Plots/"

data.matrix.LMW=read.csv("Spp_site_year matrix summarised to wetland_LMW_WL_June 2018.csv", row.names=1) # load data - object name is 'Output' ... :)

Output.LMW =data.matrix.LMW

data.dir="C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd(data.dir)
data.matrix.HLWL=read.csv("Spp_site_year matrix HTH_WL_June 2018.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.HLWL =data.matrix.HLWL


data.dir="C:/Users/jc246980/Documents/MD Vegetation/Chowilla_data_csvs/"; setwd(data.dir)
data.matrix.Chow=read.csv("Spp_site_matrix summarised to wetland_Chowilla_WL_June 2018.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.Chow =data.matrix.Chow

Output.Chow=Output.Chow[,-c(5)] # remove 'NA' column


data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
data.matrix.Gun=read.csv("Spp_site_matrix summarised to wetland_Gunbower_WL_June 2018.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.Gun =data.matrix.Gun

Output.Gun=Output.Gun[,-grep("NA",colnames(Output.Gun))] # remove 'NA' column

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/KP_data_csvs/"; setwd (data.dir) 
data.matrix.KP=read.csv("Spp_site_matrix summarised to wetland_KP_WL_June 2018.csv",row.names=1) # load data - object name is 'tada' ... :)

Output.KP =data.matrix.KP
Output.KP=Output.KP[,-grep("NA",colnames(Output.KP))] # remove 'NA' column


##########################################################################################################################################
### Step 2
### create a single matrix of species with functionalgroyups

species <- unique(c(colnames(Output.LMW),colnames(Output.HLWL),colnames(Output.Chow),colnames(Output.Gun) ,colnames(Output.KP))) # create single speies list for all sites

species=as.data.frame(sort(species))
colnames(species)="species"


# Import Hattah lakes data
data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
mydataHTWL=data.frame(read.csv("HTH_WL_FGcorrections.csv")) # functional group allocations corrected in veg database (some species had two different FG assignments)
mydataHTWL=mydataHTWL[,c("Scientific.name","Functional.group","Weed.status")]

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (data.dir)
mycodes=data.frame(read.csv("Master_Species_list.csv")) # this list has now been amended

myspecies=merge(species,mycodes,by.x="species",by.y="sp_code_simple",all.x=TRUE,all.y=FALSE)

myspecies=merge(myspecies,mydataHTWL,by.x="Scientific.name",by.y="Scientific.name",all.x=TRUE)

tada=myspecies[!duplicated(myspecies[,c('species')]),]

# Import KP data

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/KP_data_csvs/"; setwd (data.dir) 
mydataKP=data.frame(read.csv("KP_2008-2016.csv"))
mydataKPL=mydataKP[,c("Species","PFG.code.casanova","Indigenous.or.weed")]
mydataKPL=mydataKPL[!duplicated(mydataKPL[,c('Species')]),] # remove duplicates
myspecies=merge(tada,mydataKPL,by.x="Scientific.name",by.y="Species",all.x=TRUE)

my.na <- is.na(myspecies$Functional.group)
myspecies$Functional.group=as.character(myspecies$Functional.group)
myspecies$PFG.code.casanova=as.character(myspecies$PFG.code.casanova)
myspecies$Functional.group[my.na] <- myspecies$PFG.code.casanova[my.na]

myspecies$Indigenous.or.weed=as.character(myspecies$Indigenous.or.weed)
myspecies$Indigenous.or.weed <- replace(myspecies$Indigenous.or.weed, myspecies$Indigenous.or.weed=="indigenous", "FALSE")
myspecies$Indigenous.or.weed <- replace(myspecies$Indigenous.or.weed, myspecies$Indigenous.or.weed=="introduced", "TRUE")

my.na <- is.na(myspecies$Weed.status)

myspecies$Weed.status[my.na] <- myspecies$Indigenous.or.weed[my.na]

#Import LMW data


data.dir = "C:/Users/jc246980/Documents/MD Vegetation/LMW_data_csvs/"; setwd (data.dir) 
mydataLMW=data.frame(read.csv("LMW_WL.csv"))
mydataLMW=mydataLMW[,c("Scientific.name","Functional.group","Weed.status")]
mydataLMW=mydataLMW[!duplicated(mydataLMW[,c('Scientific.name')]),] # remove duplicates

myspecies1=merge(myspecies,mydataLMW,by.x="Scientific.name",by.y="Scientific.name",all.x=TRUE)

my.na <- is.na(myspecies1$Functional.group.x) # create logical vector of NAs
myspecies1$Functional.group.x=as.character(myspecies1$Functional.group.x)
myspecies1$Functional.group.y=as.character(myspecies1$Functional.group.y)
myspecies1$Functional.group.x[my.na] <- myspecies1$Functional.group.y[my.na]

my.na <- is.na(myspecies1$Weed.status.x)
myspecies1$Weed.status.x[my.na] <- myspecies1$Weed.status.y[my.na]


# Import Gun data
data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Gunbower_data_csvs/"; setwd (data.dir) 
mydataGUN=data.frame(read.csv("Gunbower_WL.csv"))
mydataGUN=mydataGUN[,c("Species","PFG","Native")]
mydataGUN=mydataGUN[!duplicated(mydataGUN[,c('Species')]),] # remove duplicates


data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (data.dir)
GunPFG=data.frame(read.csv("GunPFG.csv"))

myGUN=merge(mydataGUN,GunPFG,by.x="PFG",by.y="PFG.Code",all.x=TRUE)

myspecies2=merge(myspecies1,myGUN,by.x="Scientific.name",by.y="Species",all.x=TRUE)
my.na <- is.na(myspecies2$Functional.group.x) # create logical vector of NAs
myspecies2$PFG.y=as.character(myspecies2$PFG.y)

myspecies2$Functional.group.x[my.na] <- myspecies2$PFG.y[my.na]

my.na <- is.na(myspecies2$Weed.status.x)
myspecies2$Native=as.character(myspecies2$Native)
myspecies2$Native <- replace(myspecies2$Native, myspecies2$Native=="N", "TRUE")
myspecies2$Native <- replace(myspecies2$Native, myspecies2$Native=="Y", "FALSE")
myspecies2$Weed.status.x[my.na] <- myspecies2$Native[my.na]

write.csv(myspecies2 , file = "Species_Master_list_PFGinfo.csv") # save data out - this list has been tidied up to remove "*" from the scientifc names correct to allow a search of MCs list

###################################################################################

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (data.dir)
splist=read.csv("Species_Master_list_PFGinfo July2018.csv", row.names=1) # load data
Casanovat=read.csv("Casanova_PFGinfo.csv") # load data

myspecies3=merge(splist,Casanovat,by.x="Scientific.name.corrected",by.y="taxon",all.x=TRUE)

write.csv(myspecies3 , file = "Species_Master_list_PFG info MCPFGs.csv") # save da

####################################################################################
### Correct species names - after CC has been through and vetted PFG allocations

data.dir = "C:/Users/jc246980/Documents/MD Vegetation/Species lists for analysis/"; setwd (data.dir)
veg.df=read.csv("Species_Master_FINAL.csv") # load data

library(taxize)
library(dplyr)
library(magrittr)   

src <- c("The International Plant Names Index") # set reference 
subset(gnr_datasources(), title %in% src)# id is 167

veg.df$Scientific.name.corrected<-as.character(veg.df$Scientific.name.corrected)

result.long <- veg.df$Scientific.name.corrected %>%
gnr_resolve(data_source_ids = c(167), 
            with_canonical_ranks=T)
			
# family <-classification(veg.df$Scientific.name.corrected, db = 'itis')
# poo=family
# na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
# poo2=na.omit.list(poo)
# familyv2=as.data.frame(t(sapply(names(poo2), function (x) poo2[[x]] [,1])[c(10,11),]))



			
veg.df.APNI <- merge(veg.df,result.long,by.x="Scientific.name.corrected",by.y="user_supplied_name",all.x=TRUE)

write.csv(veg.df.APNI , file = "Species_Master_list_APNIcorrected_FINAL.csv") # save da