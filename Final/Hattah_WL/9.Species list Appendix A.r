### Appendix species list with allocations

##################################################################################################################
#### Script to generate response metrics
library(vegan)


data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)
mydata=read.csv("Spp_site_year_transect matrix HTH_WL_July 2018.csv",row.names=1) # save data out


date.dir="C:/Users/jc246980/Documents/Current projects/MD Vegetation/Environmental data/"; setwd (date.dir) 
env.data<-read.csv("Hattah Lakes wetlands transect based env data.csv") # save data out to chec

#########################################################################################################
#### Separate data into broad functional groups and exotic status

species.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Species lists for analysis/"; setwd (species.dir)
mycodes=data.frame(read.csv("Species_Master_list_APNIcorrected_BARMAHV2.csv")) # this list has now been amended
mycodes=mycodes[!duplicated(mycodes[,c('sp_code_simple','Scientific.name.corrected')]),] # remove duplicates

myspecies=as.data.frame(colnames(mydata))
colnames(myspecies)="Species"
fgrps=merge(myspecies,mycodes,by.x="Species",by.y="sp_code_simple",all.x=TRUE)

fgrps=fgrps[,c("Scientific.name.corrected","Species","Final_PFG_allocation","Weed.status")]
fgrps$Species=as.character(fgrps$Species)

data.dir = "C:/Users/jc246980/Documents/Current projects/MD Vegetation/Hattah_data_csvs/"; setwd (data.dir)

write.csv(fgrps, file = "Hattah species list Appendix A.csv") # save out for appendix