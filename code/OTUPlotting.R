# Create 3-panel figure for OTU information using underlying newt 
# data grouped by individual and labeled across disturbance
# Written by Arianna Krinos, last edits on 7 August 2018

pacman::p_load("qdapRegex", dplyr, tidyverse, reshape2)

### get the data to a usable form #####
seasonal = data.frame(read.csv("./data/LTEE_Newt_Seasonal_2.csv"))
seasonal$Phylum = unlist(rm_between(text.var = seasonal$taxonomy, left = "p__", right = ";", extract = TRUE))
seasonal$Family = unlist(rm_between(text.var = seasonal$taxonomy, left = "f__", right = ";", extract = TRUE))
seasonal$Family = gsub("\\[", "", seasonal$Family)
seasonal$Family = gsub("\\]", "", seasonal$Family)
seasonal$Family[seasonal$Family == ""] = NA

labeldata = read.csv("./data/LTEE_Newt_Seasonal_OTU_table_97_forprimer.csv")

newtdata = data.frame(read.csv("./data/newt_data_manip.csv"))
newtdata_otus = t(newtdata)
colnames(newtdata_otus) = newtdata_otus[1,]
newtdata_otus = newtdata_otus[2:nrow(newtdata_otus),]
lastrow = which(rownames(newtdata_otus) == "X") - 1

# replace with formatted data #####
newtdata_otus = seasonal
lastrow = nrow(newtdata_otus)

relativenewtdata_otus = data.frame(newtdata_otus[1:lastrow,]) %>%
  select(starts_with("LTEE"), starts_with("Phylum")) %>%
  mutate_at(vars(contains('LTEE')), as.character) %>%
  mutate_at(vars(contains('LTEE')), as.numeric) %>% replace_na(list(0)) %>%
  group_by(Phylum)  %>%
  summarize_all(funs(sum)) 
relativenewtdata_otus_melt = melt(relativenewtdata_otus)
