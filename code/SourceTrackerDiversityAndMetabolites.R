# Analyses on SourceTracker, alpha & beta diversity, and metabolite profiles, as 
# requested by Jeni on 1 October 18 via email. 
# Written by Arianna Krinos, last edits on 24 June 2020

pacman::p_load(devtools,tidyverse,stringr,qdapRegex,reshape2,plyr,RColorBrewer,vegan,proxy,philentropy)

devtools::install_github("https://github.com/GuillemSalazar/EcolUtils")
#' Effect of individual newt ID on SourceTracker results (proportion OTUs conserved)
#' 
#' Effect of individual newt ID on beta diversity (run permanova on weighted unifrac 
#' matrix with Newt ID as fixed effect)
#' 
#' Effect of individual newt ID on alpha diversity 
#' (Kruskall Wallis test on OTU richness and phylogenetic diversity)
#' 
#' Effect of individual newt ID on metabolic profiles 
#' (run permanova on Jaccard matrix with Newt ID as fixed effect)

##### Effect of individual newt ID on SourceTracker results (proportion OTUs conserved) #####
sourcetracker = read.csv("./data/LTEE_Sourcetracker_plain.csv")
newtsourcetracker = cbind(sourcetracker$AmphibID, sourcetracker$Proportion_Fall2011)
colnames(newtsourcetracker) = c("AmphibID", "Proportion_Fall2011")
kruskal.test(Proportion_Fall2011 ~ AmphibID, data = newtsourcetracker)

##### Effect of individual newt ID on beta diversity #####
unifrac = read.csv("./data/weighted_unifrac_matrix.csv")
rownames(unifrac) = unifrac[,1]
unifrac = unifrac[,-which(colnames(unifrac) == "X")]
unifrac_characteristics = matrix(0, ncol = 4, nrow = nrow(unifrac))
colnames(unifrac_characteristics) = c("NewtID","Date","Season","Disturbance")
unifrac_characteristics = as.data.frame(unifrac_characteristics)
toassess = strsplit(colnames(unifrac), '[.,"LTEE"]') # the split up for analysis pieces
timenums = c(1,2,3,4,5,6,7,8,9,11,12,13)
timepoints = c("9.9.11",
              "10.21.11",
              "12.2.11",
              "4.4.12",
              "5.16.12",
              "6.27.12",
              "8.9.12",
              "9.19.12",
              "12.5.12",
              "5.29.13",
              "8.14.13",
              "9.18.13")
disturbance = c("Pre substrate addition",
                "Pre substrate addition",
                "Pre substrate addition",
                "Pre substrate addition",
                "Pre substrate addition",
                "Pre substrate addition",
                "Post substrate addition 1",
                "Post substrate addition 1",
                "Post substrate addition 1",
                "Post substrate addition 1",
                "Post substrate addition 2",
                "Post substrate addition 2")
seasons = c("Fall","Fall","Winter","Spring","Spring","Summer","Summer","Fall","Winter","Spring","Summer","Fall")

# Go through and convert the rownames to defined characteristics
for (g in c(1:length(toassess))) {
  currentstring = toassess[[g]][toassess[[g]] != ""]
  if (currentstring[1] == "2") {
    currentstring = currentstring[-1]
  }
  amphib = as.numeric(currentstring[1])
  timepoint = which(timenums == as.numeric(currentstring[2]))
  unifrac_characteristics$NewtID[g] = amphib
  unifrac_characteristics$Date[g] = timepoints[timepoint]
  unifrac_characteristics$Season[g] = seasons[timepoint]
  unifrac_characteristics$Disturbance[g] = disturbance[timepoint]
}

adonis2(unifrac ~ NewtID, data = unifrac_characteristics)

##### Effect of disturbance on beta diversity #####
adonis2(unifrac ~ Disturbance, data = unifrac_characteristics)
pairwise.adonis(unifrac_characteristics[,c(1:3)], unifrac_characteristics$Disturbance, sim.function = "vegdist",
                sim.method = "bray")

pairwise.adonis(unifrac, unifrac_characteristics$Disturbance, sim.function = "vegdist",
                sim.method = "bray")
adonis.pair.rev <- function (dist.mat, Factor, nper = 1000, corr.method = "fdr") 
{
  require(vegan)
  as.factor(Factor)
  comb.fact <- combn(levels(Factor), 2)
  pv <- NULL
  R2 <- NULL
  for (i in 1:dim(comb.fact)[2]) {
    model.temp <- vegan::adonis2(as.dist(as.matrix(dist.mat)[Factor == comb.fact[1, i] | 
                                                               Factor == comb.fact[2, i], Factor == 
                                                       comb.fact[1, i] | 
                                                         Factor == comb.fact[2, i]]) ~ 
                                   Factor[Factor == comb.fact[1, i] | Factor == comb.fact[2, i]], 
                                 permutations = nper)
    pv <- c(pv, model.temp$`Pr(>F)`[1])
    R2 <- c(R2, model.temp$R2[1])
  }
  pv.corr <- p.adjust(pv, method = corr.method)
  data.frame(combination = paste(comb.fact[1, ], comb.fact[2,], sep = " <-> "), R2 = R2, P.value = pv, P.value.corrected = pv.corr)
}


EcolUtils::adonis.pair(vegdist(unifrac), 
                       factor(unifrac_characteristics$Disturbance))

adonis.pair.rev(vegdist(unifrac), 
            factor(unifrac_characteristics$Disturbance))

##### Effect of season on beta diversity #####
adonis2(unifrac ~ Season, data = unifrac_characteristics)
adonis.pair.rev(vegdist(unifrac), 
                factor(unifrac_characteristics$Season))

##### Effect of season on beta diversity #####
adonis2(unifrac ~ Season, data = unifrac_characteristics)
adonis.pair.rev(vegdist(unifrac), 
                factor(unifrac_characteristics$Season))

##### Effect of season on beta diversity only predisturbance #####
dist.mat = unifrac[unifrac_characteristics$Disturbance == "Pre substrate addition",
                   unifrac_characteristics$Disturbance == "Pre substrate addition"]
Factor = factor(unifrac_characteristics$Season[unifrac_characteristics$Disturbance == "Pre substrate addition"])
adonis2(unifrac[unifrac_characteristics$Disturbance == "Pre substrate addition",
                unifrac_characteristics$Disturbance == "Pre substrate addition"] ~ 
          Season[unifrac_characteristics$Disturbance == "Pre substrate addition"], 
        data = unifrac_characteristics)
adonis.pair.rev(vegdist(unifrac), 
                factor(unifrac_characteristics$Season))

##### Effect of season on beta diversity only post disturbance 1 #####
dist.mat = unifrac[unifrac_characteristics$Disturbance == "Post substrate addition 1",
                   unifrac_characteristics$Disturbance == "Post substrate addition 1"]
Factor = factor(unifrac_characteristics$Season[unifrac_characteristics$Disturbance == "Post substrate addition 1"])
adonis2(unifrac[unifrac_characteristics$Disturbance == "Post substrate addition 1",
                unifrac_characteristics$Disturbance == "Post substrate addition 1"] ~ 
          Season[unifrac_characteristics$Disturbance == "Post substrate addition 1"], 
        data = unifrac_characteristics)

##### Effect of individual newt ID on alpha diversity #####
oturichness = read.csv("./data/AlphaDiv_OTUrichness.csv")
pd = read.csv("./data/AlphaDiv_PD.csv")

kruskal.test(observed_species ~ AmphibID, data = oturichness)
kruskal.test(PD_whole_tree ~ AmphibID, data = pd)


##### Effect of individual newt ID on metabolic profiles #####
metabolicrichness = read.csv("./data/LTEE metabolite richness_2.csv")
metabolicrichness = metabolicrichness %>%
  dplyr::mutate(NewtID = as.numeric(NewtID))

simil(metabolicrichness, metabolicrichness, method = "fJaccard", pairwise = TRUE)
newmetrich = reshape2::dcast(metabolicrichness, NewtID ~ Disturbance, fun.aggregate = sum, value.var = "count")
rownames(newmetrich) = newmetrich$NewtID
newmetrich = newmetrich[,-which(colnames(newmetrich) == "NewtID")]

# build the jaccard similarity matrix using aggregated df
newestmetrich = data.frame(metabolicrichness %>%
  group_by(NewtID) %>%
  summarise_if(is.numeric, sum))

similarity_jaccard_newtid = distance(newestmetrich, method = "jaccard")
adonis2(similarity_jaccard_newtid ~ NewtID, data = newestmetrich)


# build the jaccard similarity matrix using aggregated df
newestmetrich = data.frame(metabolicrichness %>%
                             group_by(Season)) %>%
                             filter((Disturbance == "Post1")) %>%# | (Disturbance == "Post1")) %>% #"Post2") %>%
                             #summarise_if(is.numeric, sum)) %>%
                             ungroup() %>%
                             mutate(Season = factor(Season))
levels(newestmetrich$Season) <- c(1:4)
newestmetrich$Season = as.numeric(newestmetrich$Season)

similarity_jaccard_season = vegdist(newestmetrich %>% 
                                       dplyr::select(Season, count), method = "jaccard")
sampledata = read.csv("./data/sampledata.csv", row.names = 1)

adonis2(similarity_jaccard_season ~ Season, data = newestmetrich)
pairwise.adonis(newestmetrich %>% 
                  dplyr::select(count), newestmetrich$Season, sim.function = "vegdist",
                sim.method = "bray")

counts_and_season = metabolicrichness %>% select(Season,count)
levels(counts_and_season$Season) <- c(1:4)
counts_and_season$Season = as.numeric(counts_and_season$Season)

similarity_jaccard_season = distance(counts_and_season, method = "jaccard")
adonis2(similarity_jaccard_season ~ Season, data = counts_and_season)
