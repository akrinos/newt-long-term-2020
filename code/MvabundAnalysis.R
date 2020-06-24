# Script to perform and collate results easily from the mvabund analysis to
# be included in the paper 
# Written by Arianna Krinos, last edits on 8 September 2018

library(mvabund)
library(tidyverse)
library(stringr)
library(qdapRegex)
library(reshape2)
library(plyr)

seasonal = data.frame(read.csv("./data/LTEE_Newt_Seasonal.csv"))
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

numbersonly = newtdata_otus %>% select(starts_with("LTEE"))
numbersdf = mvabund(data.frame(numbersonly))
labeldata = read.csv("./data/LTEE_Newt_Seasonal_OTU_table_97_forprimer.csv")
seasons = labeldata$Season

#### Check for differences across phyla #####
# Check for patterns in our residuals
mod1 <- manyglm(numbersdf ~ newtdata_otus$Phylum, family="poisson")
plot(mod1) # fan shape present
mod1 <- manyglm(numbersdf ~ newtdata_otus$Phylum, family="negative_binomial")
plot(mod1) # no defined pattern 

anova(mod1)
anova(mod1, p.uni = "adjusted")

#### Check for differences across families #####
families <- manyglm(numbersdf ~ newtdata_otus$Family, family="poisson")
plot(families) # fan shape present
families <- manyglm(numbersdf ~ newtdata_otus$Family, family="negative_binomial")
plot(families) # no defined pattern 

anova(families)
anova(families, p.uni = "adjusted")

#### Check for differences across disturbance #####
ff = strsplit(colnames(numbersdf), "\\.")
disturbance = c()
for (t in 1:length(ff[])) {
  curr = ff[[t]]
  sample = as.numeric(curr[length(curr)])
  if (sample <= 6) {
    disturbance = c(disturbance, "Pre")
  } else if (sample <= 11) {
    disturbance = c(disturbance, "Post1")
  } else {
    disturbance = c(disturbance, "Post2")
  }
}
disturbanceanova <- manyglm(t(numbersdf) ~ disturbance, family="poisson")
plot(disturbanceanova) # fan shape present
disturbanceanova <- manyglm(t(numbersdf) ~ disturbance, family="negative_binomial")
plot(disturbanceanova) # no defined pattern 

anova(disturbanceanova)
anova(disturbanceanova, p.uni = "adjusted")


##### Check whether families are different pre-disturbance and post-disturbance 1 #####

families_grouping_tomanip = newtdata_otus %>% group_by(Family) %>%
  mutate(sums_pre = rowSums(select(numbersonly, which(disturbance == "Pre")))) %>%
  mutate(sums_post1 = rowSums(select(numbersonly, which(disturbance == "Post1")))) %>%
  mutate(sums_post2 = rowSums(select(numbersonly, which(disturbance == "Post2")))) %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(byfamily_pre = sum(sums_pre), byfamily_post1 = sum(sums_post1), 
                   byfamily_post2 = sum(sums_post2))

families_grouping_tomanip$byfamily_pre = families_grouping_tomanip$byfamily_pre / sum(families_grouping_tomanip$byfamily_pre)
families_grouping_tomanip$byfamily_post1 = families_grouping_tomanip$byfamily_post1 / sum(families_grouping_tomanip$byfamily_post1)
families_grouping_tomanip$byfamily_post2 = families_grouping_tomanip$byfamily_post2 / sum(families_grouping_tomanip$byfamily_post2)

families_grouping = data.frame(newtdata_otus %>% group_by(Family) %>%
                                 summarize_if(is.numeric, sum)) 
famnames = families_grouping$Family
families_grouping = data.frame(t(families_grouping))

colnames(families_grouping) = famnames
families_grouping = families_grouping[2:nrow(families_grouping),]

families_grouping = families_grouping %>%
  mutate(disturbance = disturbance)

famanova = families_grouping[which((families_grouping$disturbance == "Pre") | (families_grouping$disturbance == "Post1"))
                             ,1:(ncol(families_grouping)-1)]  # just pre-disturbance to post-1
famanova = apply(famanova, 2, as.numeric)

justdist1 = mvabund(famanova)
justdist1 = manyglm(famanova ~ disturbance[disturbance != "Post2"], family = "negative_binomial")
plot(justdist1)

ret = anova(justdist1, p.uni = "adjusted")

famanova = families_grouping[which((families_grouping$disturbance == "Post2") | (families_grouping$disturbance == "Post1"))
                             ,1:(ncol(families_grouping)-1)]  # just post-1 to post-2
famanova = apply(famanova, 2, as.numeric)

justdist2 = mvabund(famanova)
justdist2 = manyglm(famanova ~ disturbance[disturbance != "Pre"], family = "negative_binomial")
plot(justdist2)

ret2 = anova(justdist2, p.uni = "adjusted")

ret3 = mvabund(apply(families_grouping[,1:(ncol(families_grouping)-1)], 2, as.numeric))
ret3 = manyglm(ret3 ~ disturbance, family = "negative_binomial")
plot(ret3)

ret3 = anova(ret3, p.uni = "adjusted")


##### MVABUND FOR PHYLA #####
phyla_grouping_tomanip = newtdata_otus %>% group_by(Phylum) %>%
  mutate(sums_pre = rowSums(select(numbersonly, which(disturbance == "Pre")))) %>%
  mutate(sums_post1 = rowSums(select(numbersonly, which(disturbance == "Post1")))) %>%
  mutate(sums_post2 = rowSums(select(numbersonly, which(disturbance == "Post2")))) %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarize(byPhylum_pre = sum(sums_pre), byPhylum_post1 = sum(sums_post1), 
                   byPhylum_post2 = sum(sums_post2))

phyla_grouping_tomanip$byPhylum_pre = phyla_grouping_tomanip$byPhylum_pre / sum(phyla_grouping_tomanip$byPhylum_pre)
phyla_grouping_tomanip$byPhylum_post1 = phyla_grouping_tomanip$byPhylum_post1 / sum(phyla_grouping_tomanip$byPhylum_post1)
phyla_grouping_tomanip$byPhylum_post2 = phyla_grouping_tomanip$byPhylum_post2 / sum(phyla_grouping_tomanip$byPhylum_post2)

phyla_grouping = data.frame(newtdata_otus %>% group_by(Phylum) %>%
                              summarize_if(is.numeric, sum)) 
famnames = phyla_grouping$Phylum
phyla_grouping = data.frame(t(phyla_grouping))

colnames(phyla_grouping) = famnames
phyla_grouping = phyla_grouping[2:nrow(phyla_grouping),]

phyla_grouping = phyla_grouping %>%
  mutate(disturbance = disturbance)

famanova = phyla_grouping[which((phyla_grouping$disturbance == "Pre") | (phyla_grouping$disturbance == "Post1"))
                          ,1:(ncol(phyla_grouping)-1)]  # just pre-disturbance to post-1
famanova = apply(famanova, 2, as.numeric)

justdist1p = mvabund(famanova)
justdist1p = manyglm(famanova ~ disturbance[disturbance != "Post2"], Phylum = "negative_binomial")
plot(justdist1p)

retp = anova(justdist1p, p.uni = "adjusted")

famanova = phyla_grouping[which((phyla_grouping$disturbance == "Post2") | (phyla_grouping$disturbance == "Post1"))
                          ,1:(ncol(phyla_grouping)-1)]  # just post-1 to post-2
famanova = apply(famanova, 2, as.numeric)

justdist2p = mvabund(famanova)
justdist2p = manyglm(famanova ~ disturbance[disturbance != "Pre"], Phylum = "negative_binomial")
plot(justdist2p)

ret2p = anova(justdist2p, p.uni = "adjusted")

ret3p = mvabund(apply(phyla_grouping[,1:(ncol(phyla_grouping)-1)], 2, as.numeric))
ret3p = manyglm(ret3p ~ disturbance, Phylum = "negative_binomial")
plot(ret3p)

ret3p = anova(ret3p, p.uni = "adjusted")


#### Check for differences across Core OTUs #####
coreotu_grouping_tomanip = newtdata_otus %>% 
  mutate(sums_pre = rowSums(select(numbersonly, which(disturbance == "Pre")))) %>%
  mutate(sums_post1 = rowSums(select(numbersonly, which(disturbance == "Post1")))) %>%
  mutate(sums_post2 = rowSums(select(numbersonly, which(disturbance == "Post2")))) %>%
  dplyr::filter(X.OTU.ID %in% c("140359", "235695", "654212", "2881877", "4320368", "4378239", "4414951")) %>%
  #dplyr::select(c("X.OTU.ID", "sums_pre", "sums_post1", "sums_post2")) %>%
  dplyr::mutate(relativepre = sums_pre / sum(sums_pre), relativepost1 = sums_post1 / sum(sums_post1),
                relativepost2 = sums_post2 / sum(sums_post2))
anovaready = t(data.frame(select(coreotu_grouping_tomanip, starts_with("LTEE"))))
colnames(anovaready) = coreotu_grouping_tomanip$X.OTU.ID#c("140359", "235695", "654212", "2881877", "4320368", "4378239", "4414951")

justdist1 = mvabund(anovaready[which(disturbance == "Pre" | disturbance == "Post1"),])
firstdist = disturbance[disturbance == "Pre" | disturbance == "Post1"]
ret3p = manyglm(justdist1 ~ firstdist, family = "negative_binomial")
plot(ret3p)
coreotu = anova(ret3p, p.uni = "adjusted")

retc = manyglm(anovaready ~ disturbance, family = "negative_binomial")
plot(retc)
coreotu_c = anova(retc, p.uni = "adjusted")

justdist2 = mvabund(anovaready[which(disturbance == "Post2" | disturbance == "Post1"),])
seconddist = disturbance[disturbance == "Post2" | disturbance == "Post1"]
ret3p = manyglm(justdist2 ~ seconddist, family = "negative_binomial")
plot(ret3p)
coreotu_c2 = anova(ret3p, p.uni = "adjusted")

##### Check for differences across seasons #####
library(tidyselect)

seasonsetup = newtdata_otus %>% group_by(Family) %>%
              dplyr::summarize_if(is.numeric, funs(sums = sum))

seasondf = cbind(t(seasonsetup[,2:ncol(seasonsetup)]), seasons)
colnames(seasondf) = c(seasonsetup$Family, "Season")
seasondf = apply(data.frame(seasondf) / rowSums(seasondf), 2, as.numeric)

groupedseasonal_fam = data.frame(seasondf) %>% group_by(Season) %>%
  summarize_if(is.numeric, funs(season = sum))
percentageseasonal_fam = groupedseasonal_fam / rowSums(groupedseasonal_fam[,2:ncol(groupedseasonal_fam)])

season1 = mvabund(seasondf[,1:(ncol(seasondf) - 1)])
season1 = manyglm(season1 ~ seasondf$Season, family = "negative_binomial")
plot(season1)
season_overall = anova(season1, p.uni = "adjusted")

## check difference between season 1 and 2
seasondf_checkout = seasondf %>% filter(Season <= 2)
season2 = mvabund(seasondf_checkout[,1:(ncol(seasondf_checkout) - 1)])
season2 = manyglm(season2 ~ seasondf_checkout$Season, family = "negative_binomial")
plot(season2)
season_1_2 = anova(season2, p.uni = "adjusted")

## check difference between season 2 and 3
seasondf_checkout = seasondf %>% filter(Season <= 3 & Season >= 2)
season3 = mvabund(seasondf_checkout[,1:(ncol(seasondf_checkout) - 1)])
season3 = manyglm(season3 ~ seasondf_checkout$Season, family = "negative_binomial")
plot(season3)
season_2_3 = anova(season2, p.uni = "adjusted")

## check difference between season 3 and 4
seasondf_checkout = seasondf %>% filter(Season <= 4 & Season >= 3)
season4 = mvabund(seasondf_checkout[,1:(ncol(seasondf_checkout) - 1)])
season4 = manyglm(season4 ~ seasondf_checkout$Season, family = "negative_binomial")
plot(season4)
season_3_4 = anova(season2, p.uni = "adjusted")

## check difference between season 1 and 4
seasondf_checkout = seasondf %>% filter(Season == 4 | Season == 1)
season4 = mvabund(seasondf_checkout[,1:(ncol(seasondf_checkout) - 1)])
season4 = manyglm(season4 ~ seasondf_checkout$Season, family = "negative_binomial")
plot(season4)
season_1_4 = anova(season2, p.uni = "adjusted")

## check difference between season 2 and 4
seasondf_checkout = seasondf %>% filter(Season == 4 | Season == 2)
season4 = mvabund(seasondf_checkout[,1:(ncol(seasondf_checkout) - 1)])
season4 = manyglm(season4 ~ seasondf_checkout$Season, family = "negative_binomial")
plot(season4)
season_2_4 = anova(season2, p.uni = "adjusted")

##  By phyla this time
seasonsetup = newtdata_otus %>% group_by(Phylum) %>%
  dplyr::summarize_if(is.numeric, funs(sums = sum))

seasondf = cbind(t(seasonsetup[,2:ncol(seasonsetup)]), seasons)
colnames(seasondf) = c(seasonsetup$Phylum, "Season")
seasondf = data.frame(seasondf)

groupedseasonal = data.frame(seasondf) %>% group_by(Season) %>%
  summarize_if(is.numeric, funs(season = sum))
percentageseasonal = groupedseasonal / rowSums(groupedseasonal[,2:ncol(groupedseasonal)])s

season1 = mvabund(seasondf[,1:(ncol(seasondf) - 1)])
season1 = manyglm(season1 ~ seasondf$Season, family = "negative_binomial")
plot(season1)
season_overall = anova(season1, p.uni = "adjusted")

