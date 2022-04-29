# Groups newt data via phylum, family, and OTU and the performed mean and standard 
# deviation analysis on the available data
# Written by Arianna Krinos, last edits on 21 September 2019

pacman::p_load(tidyverse,stringr,qdapRegex,reshape2,plyr,RColorBrewer,gdata)

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

#### Just get the number for mean relative abundance of phyla #####
relativenewtdata_otus = data.frame(newtdata_otus[1:lastrow,]) %>%
                        select(starts_with("LTEE"), starts_with("Phylum")) %>%
                        mutate_at(vars(contains('LTEE')), as.character) %>%
                        mutate_at(vars(contains('LTEE')), as.numeric) %>% replace_na(list(0)) %>%
                        group_by(Phylum)  %>%
                        summarize_all(funs(sum)) 
relativenewtdata_otus_melt = melt(relativenewtdata_otus)
getrelativemeansandsds = relativenewtdata_otus_melt %>%
                             group_by(variable) %>%
                             dplyr::mutate(value = value / sum(value)) %>%
                             group_by(Phylum) %>%
                             dplyr::summarize(means = mean(value), sds = sd(value)) 
overallmeans = melt(relativenewtdata_otus) %>%
               group_by(Phylum) %>%
               dplyr::summarize(means = sum(value)) %>%
               group_by() %>%
               dplyr::mutate(means = means / sum(means))
write.csv(overallmeans, "./results/means_phyla.csv", row.names = FALSE)

#### Get the means for each of the disturbance types too for phyla #####
getdisturbancelevelresults = relativenewtdata_otus_melt %>%
                             dplyr::mutate(disturbance = as.character(newtdata$Substrate.Addition[match(as.character(relativenewtdata_otus_melt$variable), as.character(newtdata$X.OTU.ID))])) %>%
                             group_by(Phylum, disturbance) %>%
                             dplyr::summarize(means = sum(value)) %>%
                             group_by(disturbance) %>%
                             dplyr::mutate(means = means / sum(means))
getdisturbancelevelresults$disturbance =  plyr::revalue(getdisturbancelevelresults$disturbance, 
        c("Pre substrate addition"="Pre", "Post substrate addition 1"="Post-1", 
          "Post substrate addition 2" = "Post-2"))
getdisturbancelevelresults$disturbance = factor(getdisturbancelevelresults$disturbance)
phylum1 = getdisturbancelevelresults %>%
  mutate(Phylum = forcats::fct_reorder(Phylum, as.numeric(means), .fun = mean)) %>%
  ggplot(aes(x = disturbance, y = means, fill = Phylum)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ylab("Average Relative Abundance") + xlab("Disturbance") + 
  scale_x_discrete(limits=c("Pre", "Post-1", "Post-2")) + theme_classic()
write.csv(getdisturbancelevelresults, "./results/phyla_disturbance.csv", row.names = FALSE)

#### Get means for disurbance types and individuals for phyla #####
getdisturbancelevelresults = relativenewtdata_otus_melt %>%
                              dplyr::mutate(disturbance = as.character(newtdata$Substrate.Addition[match(as.character(relativenewtdata_otus_melt$variable), as.character(newtdata$X.OTU.ID))])) %>%
                              dplyr::mutate(newt = as.character(newtdata$AmphibID[match(as.character(relativenewtdata_otus_melt$variable), as.character(newtdata$X.OTU.ID))])) %>%
                              group_by(Phylum, disturbance, newt) %>%
                              dplyr::summarize(means = sum(value)) %>%
                              group_by(disturbance, newt) %>%
                              dplyr::mutate(means = means / sum(means)) %>%
                              group_by() %>%
                              dplyr::mutate(newt = factor(newt))
getdisturbancelevelresults$disturbance =  revalue(getdisturbancelevelresults$disturbance, 
                                                  c("Pre substrate addition"="Pre", "Post substrate addition 1"="Post-1", 
                                                    "Post substrate addition 2" = "Post-2"))
getdisturbancelevelresults$disturbance = factor(getdisturbancelevelresults$disturbance)
getdisturbancelevelresults$disturbance = factor(getdisturbancelevelresults$disturbance,
                                                levels(getdisturbancelevelresults$disturbance)[c(3,1,2)])
posmatches = match(sort(as.numeric(as.character(levels(getdisturbancelevelresults$newt)))), 
                   as.numeric(as.character(levels(getdisturbancelevelresults$newt))))
getdisturbancelevelresults$newt = factor(getdisturbancelevelresults$newt,
                                                levels(getdisturbancelevelresults$newt)[c(posmatches)])
phylum_plot_dist = getdisturbancelevelresults %>%
  mutate(Phylum = forcats::fct_reorder(Phylum, as.numeric(means), .fun = mean)) %>%
  ggplot(aes(x = newt, y = means, fill = Phylum)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_wrap(~disturbance, strip.position = "bottom", scales = "free_x") +
  ylab("Average Relative Abundance") + xlab("Newt") + theme_classic() + theme(legend.position = "none")
write.csv(getdisturbancelevelresults, "./results/phyla_newt_disturbance.csv", row.names = FALSE)

#### Just get the number for mean relative abundance of families #####
newtdata_otus_formax = newtdata_otus # we are looking for the most abundant of all
provis =newtdata_otus_formax %>% select(starts_with("LTEE"))
sumsRows = rowSums(provis)
newtdata_otus_formax$Family[which(sumsRows == max(sumsRows))]
sumsRows[which(sumsRows == max(sumsRows))] / sum(sumsRows)
mean(as.numeric(provis[which(sumsRows == max(sumsRows)),] / apply(provis, 2, sum)))
sd(as.numeric(provis[which(sumsRows == max(sumsRows)),] / apply(provis, 2, sum)))
relativenewtdata_otus = data.frame(newtdata_otus[1:lastrow,]) %>%
                          select(starts_with("LTEE"), starts_with("Family")) %>%
                          mutate_at(vars(contains('LTEE')), as.character) %>%
                          mutate_at(vars(contains('LTEE')), as.numeric) %>% replace_na(list(0)) %>%
                          group_by(Family)  %>%
                          summarize_all(funs(sum)) 
relativenewtdata_otus_melt = melt(relativenewtdata_otus)
getrelativemeansandsds = relativenewtdata_otus_melt %>%
                          group_by(variable) %>%
                          #mutate(sums = as.numeric(by(relativenewtdata_otus_melt$value, relativenewtdata_otus_melt$variable, sum))) %>%
                          dplyr::mutate(value = value / sum(value)) %>%
                          group_by(Family) %>%
                          dplyr::summarize(means = mean(value), sds = sd(value)) 
getrelativemeansandsds = getrelativemeansandsds[order(getrelativemeansandsds$means, decreasing = TRUE),]
getrelativemeansandsds$Family = as.character(getrelativemeansandsds$Family)
relativenewtdata_otus_melt$Family = as.character(relativenewtdata_otus_melt$Family)
relativenewtdata_otus_melt$Family[which(relativenewtdata_otus_melt$Family %in% getrelativemeansandsds$Family[19:nrow(getrelativemeansandsds)])] = "Other Families"
getrelativemeansandsds$Family[19:nrow(getrelativemeansandsds)] = "Other Families"

relativenewtdata_otus_melt$Family[relativenewtdata_otus_melt$Family == "<NA>"] = "Other Families"
relativenewtdata_otus_melt$Family[is.na(relativenewtdata_otus_melt$Family)] = "Other Families"
getrelativemeansandsds$Family[getrelativemeansandsds$Family == "<NA>" | is.na(getrelativemeansandsds$Family)] = "Other Families"
getrelativemeansandsds$means[which(getrelativemeansandsds$Family == "Other Families")[1]] = sum(getrelativemeansandsds$means[getrelativemeansandsds$Family == "Other Families"])
getrelativemeansandsds = getrelativemeansandsds[1:18,]
overallmeans = melt(relativenewtdata_otus) %>%
                group_by(Family) %>%
                dplyr::summarize(means = sum(value)) %>%
                group_by() %>%
                dplyr::mutate(means = means / sum(means))

overallmeans_otherfamilies = relativenewtdata_otus_melt %>%
                group_by(Family) %>%
                dplyr::summarize(means = sum(value)) %>%
                group_by() %>%
                dplyr::mutate(means = means / sum(means))
write.csv(overallmeans_otherfamilies, "./results/means_families.csv", row.names = FALSE)

#### Get the means for each of the disturbance types too, family #####
getdisturbancelevelresults = relativenewtdata_otus_melt %>%
                              dplyr::mutate(disturbance = as.character(newtdata$Substrate.Addition[match(as.character(relativenewtdata_otus_melt$variable), as.character(newtdata$X.OTU.ID))])) %>%
                              group_by(Family, disturbance) %>%
                              dplyr::summarize(means = sum(value)) %>%
                              group_by(disturbance) %>%
                              dplyr::mutate(means = means / sum(means))
getdisturbancelevelresults$disturbance =  revalue(getdisturbancelevelresults$disturbance, 
                                                  c("Pre substrate addition"="Pre", "Post substrate addition 1"="Post-1", 
                                                    "Post substrate addition 2" = "Post-2"))
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
barney.colors <- colorRampPalette(c("blue", "magenta", "rosybrown","darkblue", "salmon", "wheat4", "pink", "aquamarine", "yellow", "slateblue", "purple", "orange"))
split_1 = round(length(unique(getdisturbancelevelresults$Family)) / 2,0)
split_2 = length(unique(getdisturbancelevelresults$Family)) - split_1
full_rainbow = c(rainbow(split_1+split_2,rev=TRUE))
full_rainbow[2] = "purple"
full_rainbow[4] = "brown"
family_plot_dist = getdisturbancelevelresults %>%
  mutate(Family = forcats::fct_reorder(Family, as.numeric(means), .fun = mean)) %>%  
  ggplot(aes(x = disturbance, y = means, fill = Family)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ylab("Average Relative Abundance") + xlab("Disturbance") + 
  scale_fill_manual(values = c(brewer.pal(split_1,"Set3"),
                   brewer.pal(split_2,"Paired"))) + 
  scale_x_discrete(limits=c("Pre", "Post-1", "Post-2")) + theme_classic() + 
  theme(legend.key.size = unit(0.9, "line"))
write.csv(getdisturbancelevelresults, "./results/families_disturbance.csv", row.names = FALSE)

#### Get means for disurbance types and individuals, families #####
getdisturbancelevelresults = relativenewtdata_otus_melt %>%
                              dplyr::mutate(disturbance = as.character(newtdata$Substrate.Addition[match(as.character(relativenewtdata_otus_melt$variable), as.character(newtdata$X.OTU.ID))])) %>%
                              dplyr::mutate(newt = as.character(newtdata$AmphibID[match(as.character(relativenewtdata_otus_melt$variable), as.character(newtdata$X.OTU.ID))])) %>%
                              group_by(Family, disturbance, newt) %>%
                              dplyr::summarize(means = sum(value)) %>%
                              group_by(disturbance, newt) %>%
                              dplyr::mutate(means = means / sum(means)) %>%
                              group_by() %>%
                              dplyr::mutate(newt = factor(newt))
getdisturbancelevelresults$disturbance =  revalue(getdisturbancelevelresults$disturbance, 
                                                  c("Pre substrate addition"="Pre", "Post substrate addition 1"="Post-1", 
                                                    "Post substrate addition 2" = "Post-2"))
getdisturbancelevelresults$disturbance = factor(getdisturbancelevelresults$disturbance)
getdisturbancelevelresults$disturbance = factor(getdisturbancelevelresults$disturbance,
                                                levels(getdisturbancelevelresults$disturbance)[c(3,1,2)])
posmatches = match(sort(as.numeric(as.character(levels(getdisturbancelevelresults$newt)))), 
                   as.numeric(as.character(levels(getdisturbancelevelresults$newt))))
getdisturbancelevelresults$newt = factor(getdisturbancelevelresults$newt,
                                         levels(getdisturbancelevelresults$newt)[c(posmatches)])
family_plot_dist_newt =  getdisturbancelevelresults %>%
  mutate(Newt_Name = paste("Newt", match(newt,newt.ids), sep = " ")) %>%
  mutate(Family = forcats::fct_reorder(Family, as.numeric(means), .fun = mean)) %>%  
  ggplot(aes(x = disturbance, y = means, fill = Family)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_wrap(~Newt_Name, strip.position = "bottom", scales = "free_x",ncol =4) +
  scale_fill_manual(values = c(brewer.pal(split_1,"Set3"),
                               brewer.pal(split_2,"Paired"))) + 
  ylab("Average Relative Abundance") + xlab("Disturbance") + theme_classic() + theme(legend.position = "none")
family_plot_dist_newt_bydisturb =  getdisturbancelevelresults %>%
  mutate(Newt_Name = paste("Newt", match(newt,newt.ids), sep = " ")) %>%
  mutate(Family = forcats::fct_reorder(Family, as.numeric(means), .fun = mean)) %>%  
  ggplot(aes(x = newt, y = means, fill = Family)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_wrap(~disturbance, strip.position = "bottom", scales = "free_x") +
  scale_fill_manual(values = c(brewer.pal(split_1,"Set3"),
                               brewer.pal(split_2,"Paired"))) + 
  ylab("Average Relative Abundance") + xlab("Newt") + theme_classic() + theme(legend.position = "none")

write.csv(getdisturbancelevelresults, "./results/family_newt_disturbance.csv", row.names = FALSE)

# Do the seasonal-level analysis for both family and phyla #####
seasonsonly = FALSE #FALSE #TRUE
relativenewtdata_otus = data.frame(newtdata_otus[1:lastrow,]) %>%
                        select(starts_with("LTEE"), starts_with("Phylum")) %>%
                        mutate_at(vars(contains('LTEE')), as.character) %>%
                        mutate_at(vars(contains('LTEE')), as.numeric) %>% replace_na(list(0)) %>%
                        group_by(Phylum)  %>%
                        summarize_all(funs(sum)) %>%
                        melt() %>%
                        dplyr::mutate(season = labeldata$Season[match(as.character(variable), as.character(labeldata$X.OTU.ID))])%>%
                        dplyr::mutate(newt = as.character(newtdata$AmphibID[match(as.character(variable), as.character(newtdata$X.OTU.ID))])) %>%
                        group_by(Phylum, season, newt) %>%
                        dplyr::summarize(means = sum(value)) %>%
                        group_by(season, newt) %>%
                        dplyr::mutate(means = means / sum(means))
if (seasonsonly) {
  relativenewtdata_otus = data.frame(newtdata_otus[1:lastrow,]) %>%
    select(starts_with("LTEE"), starts_with("Phylum")) %>%
    mutate_at(vars(contains('LTEE')), as.character) %>%
    mutate_at(vars(contains('LTEE')), as.numeric) %>% replace_na(list(0)) %>%
    group_by(Phylum)  %>%
    summarize_all(funs(sum)) %>%
    melt() %>%
    dplyr::mutate(season = labeldata$Season[match(as.character(variable), as.character(labeldata$X.OTU.ID))])%>%
    group_by(Phylum, season) %>%
    dplyr::summarize(means = sum(value)) %>%
    group_by(season) %>%
    dplyr::mutate(means = means / sum(means))
  
  relativenewtdata_otus$season = factor(relativenewtdata_otus$season)
  relativenewtdata_otus$season = factor(relativenewtdata_otus$season,
                                                  levels(relativenewtdata_otus$season)[c(2,3,1,4)])
  phylum_plot_season = relativenewtdata_otus %>%
    mutate(Phylum = forcats::fct_reorder(Phylum, as.numeric(means), .fun = mean)) %>%
    ggplot(aes(x = season, y = means, fill = Phylum)) + 
    geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("Average Relative Abundance") + xlab("Season") + theme_classic()
  write.csv(relativenewtdata_otus, "./results/phyla_season_newt.csv", row.names = FALSE)
} else {
    relativenewtdata_otus$newt = factor(relativenewtdata_otus$newt)
    posmatches = match(sort(as.numeric(as.character(levels(relativenewtdata_otus$newt)))), 
                       as.numeric(as.character(levels(relativenewtdata_otus$newt))))
    relativenewtdata_otus$newt = factor(relativenewtdata_otus$newt,
                                             levels(relativenewtdata_otus$newt)[c(posmatches)])
    
    relativenewtdata_otus$season = factor(relativenewtdata_otus$season)
    relativenewtdata_otus$season = factor(relativenewtdata_otus$season,
                                          levels(relativenewtdata_otus$season)[c(2,3,1,4)])
    phylum_plot_newt_season = relativenewtdata_otus %>%
      mutate(Phylum = forcats::fct_reorder(Phylum, as.numeric(means), .fun = mean)) %>%  
      ggplot(aes(x = newt, y = means, fill = Phylum)) + 
      geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      facet_wrap(~season, strip.position = "bottom", scales = "free_x") +
      ylab("Average Relative Abundance") + xlab("Season") + theme_classic() + theme(legend.position = "none")
    write.csv(relativenewtdata_otus, "")
}

pacman::p_load(gridExtra, cowplot)

library(cowplot)
p <- plot_grid(phylum1, phylum_plot_dist, family_plot_dist, family_plot_dist_newt, nrow = 2, labels=c('A', 'B', 'C', 'D'))
grid.arrange(phylum1, phylum_plot_dist, family_plot_dist, family_plot_dist_newt, nrow = 2)


p <- plot_grid(family_plot_dist_newt+theme(text=element_text(size=16))+xlab(""), 
               family_plot_dist+theme(text=element_text(size=16)), nrow = 2,
               labels=c('A', 'B'),rel_heights = c(0.58,0.5))
ggsave(plot=p,file="figure3.tiff",dpi=500,width=11,height=11,units="in")

seasonsonly = FALSE
relativenewtdata_otus = data.frame(newtdata_otus[1:lastrow,]) %>%
  select(starts_with("LTEE"), starts_with("Family")) %>%
  mutate_at(vars(contains('LTEE')), as.character) %>%
  mutate_at(vars(contains('LTEE')), as.numeric) %>% replace_na(list(0)) %>%
  group_by(Family)  %>%
  summarize_all(funs(sum)) %>%
  melt() %>%
  dplyr::mutate(season = as.character(labeldata$Season[match(as.character(variable), as.character(labeldata$X.OTU.ID))]))%>%
  dplyr::mutate(newt = as.character(labeldata$AmphibID[match(as.character(variable), as.character(labeldata$X.OTU.ID))]))%>%
  group_by(Family, season, newt) %>%
  dplyr::summarize(means = sum(value)) %>%
  group_by(season, newt) %>%
  dplyr::mutate(means = means / sum(means))
if (seasonsonly) {
  relativenewtdata_otus = data.frame(newtdata_otus[1:lastrow,]) %>%
    select(starts_with("LTEE"), starts_with("Family")) %>%
    mutate_at(vars(contains('LTEE')), as.character) %>%
    mutate_at(vars(contains('LTEE')), as.numeric) %>% replace_na(list(0)) %>%
    group_by(Family)  %>%
    summarize_all(funs(sum)) %>%
    melt() %>%
    dplyr::mutate(season = labeldata$Season[match(as.character(variable), as.character(labeldata$X.OTU.ID))])%>%
    group_by(Family, season) %>%
    dplyr::summarize(means = sum(value)) %>%
    group_by(season) %>%
    dplyr::mutate(means = means / sum(means))
  
  relativenewtdata_otus$season = factor(relativenewtdata_otus$season)
  relativenewtdata_otus$season = factor(relativenewtdata_otus$season,
                                        levels(relativenewtdata_otus$season)[c(2,3,1,4)])
  relativenewtdata_otus = relativenewtdata_otus %>%
    mutate(Family = forcats::fct_reorder(Family, as.numeric(means), .fun = mean)) 
  ordered = as.character(na.exclude(unique(relativenewtdata_otus$Family[order(relativenewtdata_otus$means, decreasing = TRUE)])))
  exclude = ordered[18:length(ordered)]
  notexclude = ordered[-match(exclude, ordered)]
    
  rrrr = relativenewtdata_otus %>%  
    group_by(season) %>%
    filter(Family %in% exclude | is.na(Family)) %>%
    group_by(season) %>%
    dplyr::summarize(means = sum(means)) %>%
    mutate(Family = "Other Families")
  family_plot_season =  relativenewtdata_otus %>%
    bind_rows(rrrr) %>%
    ungroup() %>%
    filter(Family %in% notexclude | Family == "Other Families") %>%
    mutate(season = factor(season)) %>%
    mutate(Family = forcats::fct_reorder(Family, as.numeric(means), .fun = mean)) %>%
    ggplot(aes(x = season, y = means, fill = Family)) + 
    scale_fill_manual(values = c(brewer.pal(split_1,"Set3"),
                                 brewer.pal(split_2,"Paired"))) + 
    geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("Average Relative Abundance") + xlab("Season") + theme_classic() + theme(legend.key.size = unit(0.9, "line"))
} else {
  relativenewtdata_otus$newt = factor(relativenewtdata_otus$newt)
  posmatches = match(sort(as.numeric(as.character(levels(relativenewtdata_otus$newt)))), 
                     as.numeric(as.character(levels(relativenewtdata_otus$newt))))
  relativenewtdata_otus$newt = factor(relativenewtdata_otus$newt,
                                      levels(relativenewtdata_otus$newt)[c(posmatches)])
  
  relativenewtdata_otus$season = factor(relativenewtdata_otus$season)
  relativenewtdata_otus$season = factor(relativenewtdata_otus$season,
                                        levels(relativenewtdata_otus$season)[c(2,3,1,4)])
  relativenewtdata_otus = relativenewtdata_otus %>%
    mutate(Family = forcats::fct_reorder(Family, as.numeric(means), .fun = mean)) 
  
  ordered = as.character(na.exclude(unique(relativenewtdata_otus$Family[order(relativenewtdata_otus$means, decreasing = TRUE)])))
  exclude = ordered[18:length(ordered)]
  notexclude = ordered[-match(exclude, ordered)]
    
  rrrr = relativenewtdata_otus %>%  
    group_by(season,newt) %>%
    filter(Family %in% exclude | is.na(Family)) %>%
    group_by(season,newt) %>%
    dplyr::summarize(means = sum(means)) %>%
    mutate(Family = "Other Families") 
  
  newt.ids = c(5,6,7,8,9,13,16,17)
  family_plot_newt_season = relativenewtdata_otus %>%
    bind_rows(rrrr) %>%
    ungroup() %>%
    filter(Family %in% notexclude | Family == "Other Families") %>%
    mutate(season = factor(season)) %>%
    #mutate(season = factor(season,
    #                       levels(relativenewtdata_otus$season)[c(2,3,1,4)])) %>%
    mutate(Family = forcats::fct_reorder(Family, as.numeric(means), .fun = mean)) %>%
    mutate(Newt_Name = paste("Newt", match(newt,newt.ids), sep = " "))%>%
    ggplot(aes(x = season, y = means, fill = Family)) + 
    geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_fill_manual(values = c(brewer.pal(split_1,"Set3"),
                                 brewer.pal(split_2,"Paired"))) + 
    facet_wrap(~Newt_Name, strip.position = "bottom", scales = "free_x", ncol = 4) +
    ylab("Average Relative Abundance") + xlab("Season") + theme_classic() + theme(legend.position = "none") + 
    theme(legend.key.size = unit(0.9, "line"))
}

plot_grid(phylum_plot_season, phylum_plot_newt_season, family_plot_season, family_plot_newt_season, 
             nrow = 2, labels = c("A", "B", "C", "D"))

plot_grid(family_plot_newt_season + xlab("")+theme(text = element_text(size=16)), 
          family_plot_season+theme(text = element_text(size=16)),
          nrow = 2, labels = c("A", "B"), rel_heights = c(0.58,0.5))
ggsave(file="figure4.tiff",dpi=500,width=11,height=11,units="in")

##### Plot simpson evenness, OTU richness, and proportion of OTUs shared with initial Fall 2011 samples #####

newtda = data.frame(read.csv("./data/newt_data_manip.csv"))
counts = apply(data.frame(apply(newtda[1:83,2:225], 2, function(x) as.numeric(as.character(x))) > 0), 1, sum)
newtda = newtda %>%
         mutate_at(vars(contains('LTEE')), as.numeric) %>% replace_na(list(0)) %>%
         mutate(OTUrich = c(counts, rep(0, nrow(newtda) - length(counts)))) %>%
         select(Date, AmphibID, OTUrich) %>%
         mutate(Date = as.Date(Date, format = "%m.%d.%y")) %>%
         filter(OTUrich != 0 & !is.na(OTUrich)) %>%
         group_by(AmphibID)

meltednewtda = melt(newtda, id.vars = c("AmphibID", "Date"))
meanmelted = meltednewtda %>%
             dplyr::group_by(Date) %>%
             dplyr::summarize(means = mean(value), sds = sd(value))

meltednewtda = meltednewtda %>%
              full_join(meanmelted)

#### OTU Richness plot #####
pacman::p_load(RColorBrewer,lubridate)
start_dates_summer = c(meltednewtda$Date[1], as.Date("6.21.12", format = "%m.%d.%y") + c(0:1) * 365)
start_dates_fall = c(as.Date("9.21.11", format = "%m.%d.%y") + c(0:1) * 365) # ends at 9/18/2013
start_dates_winter = c(as.Date("12.21.11", format = "%m.%d.%y") + c(0:1) * 365) # ends at 9/18/2013
start_dates_spring = c(as.Date("3.21.12", format = "%m.%d.%y") + c(0:1) * 365) # starts past 2011 spring

end_dates_summer = c(as.Date("9.20.11", format = "%m.%d.%y") + c(0:1) * 365, max(as.Date(meltednewtda$Date)))
end_dates_fall = c(as.Date("12.20.11", format = "%m.%d.%y") + c(0:1) * 365)
end_dates_winter = c(as.Date("3.20.12", format = "%m.%d.%y") + c(0:1) * 365)
end_dates_spring = c(as.Date("6.20.12", format = "%m.%d.%y") + c(0:1) * 365)
tomapto = unique(sort(meltednewtda$AmphibID))
for (curr in 1:nrow(meltednewtda)) {
  meltednewtda$AmphibID[curr] <- which(meltednewtda$AmphibID[curr] == tomapto)
}
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
gg = ggplot(meltednewtda) + 
  scale_x_date(date_labels="%b %Y",date_breaks  ="1 month") +
  ylab("OTU Richness") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_rect(data=data.frame(xmin=start_dates_spring,
                   xmax=end_dates_spring,
                   ymin=-Inf,
                   ymax=Inf),
   aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Spring"),
   alpha=0.3) +
  geom_rect(data=data.frame(xmin=start_dates_summer,
                            xmax=end_dates_summer,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Summer"),
            alpha=0.3) +
  geom_rect(data=data.frame(xmin=start_dates_fall,
                            xmax=end_dates_fall,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Fall"),
            alpha=0.3) +
  geom_rect(data=data.frame(xmin=start_dates_winter,
                            xmax=end_dates_winter,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Winter"),
            alpha=0.3) +
  geom_line(data = meanmelted, aes(x = (Date), y = means, group = 1), color = "black", size = 1.5) + ylim(0, 175) + 
  geom_errorbar(data = meanmelted, aes(x=(Date), ymin=means-sds, ymax=means+sds), colour="black", width=40, size = 0.75) + 
  geom_ribbon(aes(x=Date, y=means, ymax=means+sds, ymin=means-sds), 
              alpha=0.4, fill = "black") + 
  
  geom_vline(xintercept = as.Date("6.28.12", format = "%m.%d.%y"), linetype = "dotdash", size = 1.5) + 
  geom_vline(xintercept = as.Date("6.1.13", format = "%m.%d.%y"), linetype = "dotdash", size = 1.5) +
  geom_point(data = meanmelted, aes(x = (Date), y = means, color = "Mean"), size = 3.5) + 
  #geom_line(aes(x = Date, y = value, group = AmphibID, colour = factor(AmphibID)), size = 0.5) +
  geom_point(aes(x = (Date), y = value, group = AmphibID, colour = factor(AmphibID)), size = 2) +
  scale_color_manual("Newt ID", #breaks = c(1:8,"Mean"), values =c(brewer.pal(8, "Spectral"), "black")) +
                     breaks = c(levels(factor(meltednewtda$AmphibID)), "Mean"), 
                     values = safe_colorblind_palette[1:9])+
                     #values =c(brewer.pal(8, "Dark2"), "black"))+
  xlab("Date") + 
  scale_fill_manual('Season',
                    values = gray.colors(4, start = 0.2, end = 0.9),
                    breaks = c("Spring","Summer","Fall","Winter"),
                    #values = c('green', 'red', 'brown', 'lightblue'),  
                    guide = guide_legend(override.aes = list(alpha = 0.2))) 
gg
ggsave("testnewtplot.png", dpi=5000, width=8,height=6,units="in")

newtda = data.frame(read.csv("./data/newt_data_manip.csv"))
newtda = newtda %>%
  dplyr::mutate_at(vars(contains('LTEE')), as.numeric) %>% replace_na(list(0)) %>%
  dplyr::select(Date, AmphibID, simpson) %>%
  dplyr::mutate(Date = as.Date(Date, format = "%m.%d.%y")) %>%
  dplyr::filter(simpson != 0 & !is.na(simpson)) %>%
  dplyr::group_by(AmphibID)

meltednewtda = melt(newtda, id.vars = c("AmphibID", "Date"))
meanmelted = meltednewtda %>%
  dplyr::group_by(Date) %>%
  dplyr::summarize(means = mean(value), sds = sd(value)) 

##### Simpson Evenness Plot #####
hh = ggplot(meltednewtda) + 
  scale_x_date(date_labels="%b %Y",date_breaks  ="1 month") +
  ylab("Simpson Evenness") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_rect(data=data.frame(xmin=start_dates_spring,
                            xmax=end_dates_spring,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Spring"),
            alpha=0.3) +
  geom_rect(data=data.frame(xmin=start_dates_summer,
                            xmax=end_dates_summer,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Summer"),
            alpha=0.3) +
  geom_rect(data=data.frame(xmin=start_dates_fall,
                            xmax=end_dates_fall,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Fall"),
            alpha=0.3) +
  geom_rect(data=data.frame(xmin=start_dates_winter,
                            xmax=end_dates_winter,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Winter"),
            alpha=0.3) + 
  geom_ribbon(data = meanmelted, aes(x=Date, y=means, ymax=means+sds, ymin=means-sds), 
              alpha=0.4, fill = "black") + 
  geom_vline(xintercept = as.Date("6.28.12", format = "%m.%d.%y"), linetype = "dotdash", size = 1.5) + 
  geom_vline(xintercept = as.Date("6.1.13", format = "%m.%d.%y"), linetype = "dotdash", size = 1.5) +
  geom_line(data = meanmelted, aes(x = Date, y = means, group = 1), color = "black", size = 1.5) + ylim(0, 1) +
  geom_errorbar(data = meanmelted, aes(x=Date, ymin=means-sds, ymax=means+sds), colour="black", width=40, size = 1) + 
  geom_point(data = meanmelted, aes(x = Date, y = means, color = "Mean"), size = 1.5) + 
  #geom_line(aes(x = Date, y = value, group = AmphibID, colour = factor(AmphibID)), size = 1) +
  geom_point(aes(x = Date, y = value, group = AmphibID, colour = factor(AmphibID)), size = 2) +
  scale_color_manual("Newt ID", breaks = c(levels(factor(meltednewtda$AmphibID)), "Mean"),
                     values = safe_colorblind_palette[1:9])+
                     
                     #values =c(brewer.pal(8, "Dark2"), "black"))+
  scale_fill_manual('Season',
                    values = gray.colors(4, start = 0.2, end = 0.9),
                    breaks = c("Spring","Summer","Fall","Winter"),
                    #values = c('green', 'red', 'brown', 'lightblue'),  
                    guide = guide_legend(override.aes = list(alpha = 0.2))) 
hh


#### Calculating OTUs shared with original ####
newtda = data.frame(read.csv("./data/newt_data_manip.csv"))
newtda_orig = data.frame(read.csv("./data/newt_data_manip.csv"))

numericdata = apply(newtda[1:83,2:225], 2, function(x) as.numeric(as.character(x)))
newtda = aggregate(numericdata, by = newtda[1:83, c("Year", "Season", "AmphibID")], sum, na.rm = TRUE)

otusdata = t(data.frame(apply(newtda[,4:ncol(newtda)], 2, function(x) as.numeric(as.character(x))) > 0))

fullnames = paste0(newtda$Season, "_", newtda$Year, "_", newtda$AmphibID)
colnames(otusdata) = paste0(newtda$Season, "_", newtda$Year, "_", newtda$AmphibID)
otusdata = otusdata[,grep("Fall_2011_", fullnames)]
amphibids = newtda$AmphibID[grep("Fall_2011_", fullnames)]
otupresence = matrix(0, nrow = length(unique(amphibids)), ncol = nrow(otusdata))
colnames(otupresence) = rownames(otusdata)
rownames(otupresence) = sort(unique(amphibids), decreasing = FALSE)
otupresence[t(otusdata) == TRUE] = 1

counts = c(rep(0, nrow(newtda) - 6)); seasons = c("Spring", "Summer", "Fall", "Winter")
for (g in 1:length(counts)) {
  currow = colnames(newtda)[which(as.numeric(as.character(unlist(newtda[g,4:ncol(newtda)], use.names = FALSE))) > 0) + 3]
  selectedrows = colnames(otupresence)[otupresence[which(rownames(otupresence) == newtda$AmphibID[g]),] == 1]
  correct = 0; SeasonCheck = "Winter"; YearCheck = 2011; countup = 1
  while (length(selectedrows) == 0) {
    allotusdata = t(data.frame(apply(newtda[,4:ncol(newtda)], 2, function(x) as.numeric(as.character(x))) > 0))
    colnames(allotusdata) = paste0(newtda$Season, "_", newtda$Year, "_", newtda$AmphibID)
    otusdata = allotusdata[,grep(paste0(SeasonCheck,"_", as.character(YearCheck)), fullnames)]
    selectedrows = rownames(otusdata)[which(otusdata[,grep(paste0("_", newtda$AmphibID[g]), colnames(otusdata))])]
    if (countup %% 4 == 1) {
      YearCheck = YearCheck + 1
    } else {
      SeasonCheck = seasons[(which(seasons == SeasonCheck) + 1) %% 4]
    }
    countup = countup + 1
  }
  for (t in currow) {
    if (t %in% selectedrows) {
      correct = correct + 1
    }
  }
  counts[g] = correct / length(selectedrows)
}
conewtda = newtda %>%
  mutate(PercentOTU = c(counts, rep(0, nrow(newtda) - length(counts)))) %>%
  filter(!is.na(PercentOTU)) %>%
  mutate(YearSeason = paste0(Year, Season)) %>%
  select(YearSeason, AmphibID, PercentOTU) %>%
  filter(PercentOTU != 0) %>%
  group_by(AmphibID)

meltednewtda = melt(conewtda, id.vars = c("AmphibID", "YearSeason"))
meanmelted = meltednewtda %>%
  dplyr::group_by(YearSeason) %>%
  dplyr::summarize(means = mean(value), sds = sd(value))

meltednewtda = meltednewtda %>%
  full_join(meanmelted)

ii = ggplot(meltednewtda) + geom_line(aes(x = YearSeason, y = value, group = AmphibID, colour = factor(AmphibID)), size = 2) +
  geom_point(aes(x = YearSeason, y = value, group = AmphibID, colour = factor(AmphibID)), size = 2) +
  #geom_vline(xintercept = as.Date("6.28.12", format = "%m.%d.%y"), linetype = "dotted", size = 2) + 
  #geom_vline(xintercept = as.Date("6.1.13", format = "%m.%d.%y"), linetype = "dotted", size = 2) +
  geom_line(data = meanmelted, aes(x = YearSeason, y = means, group = 1), color = "black", size = 1.5) + ylim(0, 1) + 
  geom_ribbon(data = meanmelted, aes(x=YearSeason, y=means, ymax=means+sds, ymin=means-sds), 
              alpha=0.4, fill = "black") + 
  geom_errorbar(data = meanmelted, aes(x=YearSeason, ymin=means-sds, ymax=means+sds), colour="black", width=0.5, size = 1.5) + 
  geom_point(data = meanmelted, aes(x = YearSeason, y = means, color = "Mean"), size =1.5) + 
  scale_color_manual("Newt ID", breaks = c(levels(factor(meltednewtda$AmphibID)), "Mean"), 
                     values = safe_colorblind_palette[1:9])+
                     #values =c(brewer.pal(8, "Spectral"), "black"))+
  ylab("Proportion conserved Fall 2011 OTUs") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ii

# Same information but not pooled via Jeni's sourcetracker code
sourcetrack = read.csv("./data/LTEE_Sourcetracker_plain.csv")
sourcetrack = sourcetrack %>%
              dplyr::select(AmphibID, Season, Year, SeasonYear, Date, Proportion_Fall2011, Proportion_SD_Fall2011) %>%
              dplyr::filter(!is.na(AmphibID))

meanmelted = sourcetrack %>%
             dplyr::group_by(Date) %>%
             dplyr::summarize(means = mean(na.exclude(Proportion_Fall2011)), sds = sd(na.exclude(Proportion_Fall2011))) %>%
             dplyr::mutate(Date = as.Date(Date, format = "%m.%d.%y"))
orderedmm = meanmelted[order(meanmelted$Date),] # order the means by date 

ii = ggplot(sourcetrack) + 
  ylab("Proportion conserved Fall 2011 OTUs") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
  scale_x_date(date_labels="%b %Y",date_breaks  ="1 month") +
  geom_rect(data=data.frame(xmin=start_dates_spring,
                            xmax=end_dates_spring,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Spring"),
            alpha=0.3) +
  geom_rect(data=data.frame(xmin=start_dates_summer,
                            xmax=end_dates_summer,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Summer"),
            alpha=0.3) +
  geom_rect(data=data.frame(xmin=start_dates_fall,
                            xmax=end_dates_fall,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Fall"),
            alpha=0.3) +
  geom_rect(data=data.frame(xmin=start_dates_winter,
                            xmax=end_dates_winter,
                            ymin=-Inf,
                            ymax=Inf),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill = "Winter"),
            alpha=0.3) +
  geom_ribbon(data = meanmelted, aes(x=Date, y=means, ymax=means+sds, ymin=means-sds), 
              alpha=0.4, fill = "black") + 
  geom_line(data = meanmelted, aes(x = Date, y = means, group = 1), color = "black", size = 1.5) + ylim(0, 1.25) + 
  geom_vline(xintercept = as.Date("6.28.12", format = "%m.%d.%y"), linetype = "dotdash", size = 1.5) + 
  geom_vline(xintercept = as.Date("6.1.13", format = "%m.%d.%y"), linetype = "dotdash", size = 1.5) +
  geom_errorbar(data = orderedmm, aes(x=Date, ymin=means-sds, ymax=means+sds), colour="black", width=40, size = 1) + 
  geom_point(data = meanmelted, aes(x = Date, y = means, color = "Mean"), size = 1.5) +
  #geom_line(aes(x = as.Date(Date, format = "%m.%d.%y"),y = Proportion_Fall2011, group = AmphibID, colour = factor(AmphibID)), size = 1) + 
  geom_point(aes(x = as.Date(Date, format = "%m.%d.%y"), y = Proportion_Fall2011, group = AmphibID, colour = factor(AmphibID)), size = 2) +
  scale_color_manual("Newt ID", breaks = c(levels(factor(meltednewtda$AmphibID)), "Mean"), 
                     
                     values = safe_colorblind_palette[1:9])+
                     #values =c(brewer.pal(8, "Dark2"), "black"))+
  scale_fill_manual('Season',
                    values = gray.colors(4, start = 0.2, end = 0.9),
                    breaks = c("Spring","Summer","Fall","Winter"),
                    #values = c('green', 'red', 'brown', 'lightblue'),  
                    guide = guide_legend(override.aes = list(alpha = 0.2))) +
  xlab("Date")
ii

gg
hh
ii
pacman::p_load(cowplot,gridExtra,ggpubr)
legend = get_legend(gg)
moredec <- function(x) sprintf("%.2f", x)
chart = plot_grid(gg + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()),
          hh + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()),
          ii + theme(legend.position = "none") + ylab("Ratio Fall 2011 OTUs conserved") + scale_y_continuous(labels=moredec),
          nrow = 3, rel_heights = c(1,1,1.5), labels = c("A","B","C"))
final = plot_grid(chart,
          legend, rel_widths = c(2, .25), nrow = 1)
grid.arrange(gg,hh,ii,nrow = 3)
ggarrange(gg,hh,ii,nrow = 3, common.legend = TRUE, legend = "right")
ggarrange(ii,gg,hh,nrow = 3, common.legend = TRUE, legend = "right")

## Reordering on 22 Aug as per Jeni's recommendation

chart = plot_grid(ii + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()) + ylab("Proportion of Initial OTUs"),
                  gg + theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank()),
                  hh + theme(legend.position = "none"), #+ ylab("Ratio Fall 2011 OTUs conserved") + scale_y_continuous(labels=moredec),
                  nrow = 3, rel_heights = c(1,1,1.3), labels = c("A","B","C"))

final = plot_grid(chart,
                  legend, rel_widths = c(2, .25), nrow = 1)
ggsave(file="figure1.tiff",dpi=500,width=8,height=8,units="in")

ggsave("plots/Figure1.pdf",width=8.5,height=11,units="in",dpi=1600)

pdf(file = "plots/Figure1.pdf", width = 8.5, height = 11)
dev.off()

