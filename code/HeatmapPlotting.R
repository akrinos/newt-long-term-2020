# Simple code to make core OTU heatmap
# Written by Arianna Krinos, last edits on 24 June 2020

pacman::p_load(tidyverse,stringr,qdapRegex,reshape2,plyr,RColorBrewer)

newtdata_otus_cols = read.csv("./data/LTEE_Newt_Seasonal_OTU_table_97_forprimer.csv")
coreOTUs = c(140359,235695,654212,2881877,4320368,4378239,4414951)
name_coreOTUs = paste0("X",coreOTUs)
newt.ids = c(5,6,7,8,9,13,16,17)
formalOTUnames = c("Comamonadaceae 140359",
                   "Cellulomonas 235695",
                   "Methylophilaceae 654212",
                   "Cellulomonadaceae 2881877",
                   "Comamonadaceae 4320368",
                   "Hydrogenophaga 4378239",
                   "Comamonadaceae 4414951")

newtdata_otus_cols = newtdata_otus_cols %>%
  dplyr::select(c(name_coreOTUs,Substrate.Addition,Season,AmphibID,Date)) %>%
  mutate(Newt_Name = match(AmphibID,newt.ids))
newtdata_otus_cols_melt = melt(newtdata_otus_cols,
                               id.vars = c("Date","Newt_Name","Substrate.Addition",
                                           "AmphibID","Season"),
                               value.name = "Abundance") %>%
  dplyr::rename("CoreOTU"=variable,"Newt.Name"=Newt_Name) %>%
  dplyr::mutate(CoreOTU = formalOTUnames[match(CoreOTU, name_coreOTUs)])

date_Frame = cbind(unique(newtdata_otus_cols_melt$Date), 
      sort(unique(as.Date(newtdata_otus_cols_melt$Date, format="%m.%d.%y"))),
      1:length(unique(as.Date(newtdata_otus_cols_melt$Date, format="%m.%d.%y"))))
colnames(date_Frame) = c("Date","SortedDate", "SampleDate")
modded_newtdata = newtdata_otus_cols_melt %>% dplyr::left_join(date_Frame, by = "Date", copy=TRUE)
ggplot(modded_newtdata,aes(x = as.numeric(SampleDate),#as.Date(Date,format = "%m.%d.%y"),
                                   y=factor(Newt.Name),
                                   fill =as.numeric(as.character(Abundance))/29000))+
  geom_tile()+
  facet_wrap(~CoreOTU,ncol = 2) + 
  #geom_vline(aes(xintercept=7),color = "black",size =1.5)+
  #geom_vline(aes(xintercept=as.Date("8.9.12",format="%m.%d.%y")),color = "black",size =1.5)+
  geom_vline(aes(xintercept=6.5),color = "black",size =1.5,linetype="dotdash")+
  geom_vline(aes(xintercept=10.5),color = "black",size =1.5,linetype="dotdash")+
  scale_fill_gradient(low = "blue", high = "red", name = "Abundance") + 
  scale_x_continuous(labels = as.Date(as.numeric(unique(modded_newtdata$SortedDate)),origin="1970-01-01"), 
                     breaks = as.numeric(unique(modded_newtdata$SampleDate)))+#, breaks = modded_newtdata$SampleDate) + 
  ylab("Newt") + xlab("Date") + theme_classic()+
  theme(axis.text.x = element_text(angle = 75, hjust = 1),
        text = element_text(size=16))

ggsave(file="figure5.tiff",dpi=500,width=10,height=12,units="in")

  