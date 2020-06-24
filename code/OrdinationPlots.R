# Simple code to create ordination plots using the vegan package
# Written by Arianna Krinos, last edits on 24 June 2020

pacman::p_load(phyloseq,dplyr,ggplot2,vegan,adonis,ape,philentropy,loop,MASS)

##### With the full dataset - applies to disturbance and newt ##### 
otumat_newt = read.csv("./data/otutable.csv",row.names = 1)
taxmat_newt = as.matrix(read.csv("./data/taxtable.csv",row.names = 1))

OTU = otu_table(otumat_newt[,1:(ncol(otumat_newt)-1)], taxa_are_rows = TRUE)
TAX = tax_table(taxmat_newt)

physeq = phyloseq(OTU, TAX)

sampledata = read.csv("./data/sampledata.csv", row.names = 1)

random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

physeq2 = phyloseq(OTU, TAX, sampledata[match(rownames(sampledata), colnames(otumat_newt)[1:(ncol(otumat_newt)-1)]),], random_tree)

unifrac_newt = UniFrac(physeq2, weighted=TRUE)
unifrac_newt = dist(t(otumat_newt[,1:(ncol(otumat_newt)-1)]), method = "maximum")
unifrac_newt = data.frame(read.csv("data/weightedunifrac.csv",row.names = 1))

newt_NMDS=metaMDS(unifrac_newt,k=2,trymax=10000)

# from https://chrischizinski.github.io/rstats/vegan-ggplot2/
data.scores <- as.data.frame(scores(newt_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- as.character(sampledata$Disturbance[match(rownames(data.scores),rownames(sampledata))])  #  add the grp variable created earlier
head(data.scores)  #look at the data

disturb = ggplot() +  
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  theme(legend.position = "none") + 
  ylim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  xlim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  scale_colour_manual(name = "", values=c("Pre-disturbance" = "red", "Post disturbance 1" = "blue",
                               "Post disturbance 2" = "green"),
                      breaks = c("Pre-disturbance","Post disturbance 1",
                                 "Post disturbance 2"),
                      labels = c("Pre-disturbance","Post disturbance 1",
                                 "Post disturbance 2")) +
  scale_shape_manual(name = "", values = c("Pre-disturbance" = 15,
                                "Post disturbance 1" = 16,
                                "Post disturbance 2" = 17),
                     breaks = c("Pre-disturbance","Post disturbance 1",
                                "Post disturbance 2"))+
  annotate("text", x = max(data.scores$NMDS1,data.scores$NMDS2) - 0.07, 
           y = min(data.scores$NMDS2,data.scores$NMDS1) + 0.01,
             #max(data.scores$NMDS1) - 0.07, y = max(data.scores$NMDS2) - 0.01, 
           label = paste0("Stress=",as.character(round(newt_NMDS$stress,3)))) +
  theme_bw() + ylab("") + xlab("")+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.position = c(0.10,0.01),#c(0.97, 0.99),
        legend.justification = c("left", "bottom"),legend.background = element_rect(fill = NA)) #,color="black"))
  coord_equal()

##### by season for pre disturb##### 
otumat_newt = read.csv("./data/otutable.csv",row.names = 1)
taxmat_newt = as.matrix(read.csv("./data/taxtable.csv",row.names = 1))
sampledata = read.csv("./data/sampledata.csv")#, row.names = 1) 
sampledata = sampledata %>% dplyr::filter(Disturbance == "Pre-disturbance")
sampledata$X.OTU.ID = as.character(sampledata$X.OTU.ID)
otumat_newt = otumat_newt[,as.numeric(as.character(na.exclude(match(sampledata$X.OTU.ID,colnames(otumat_newt)))))]

OTU = otu_table(otumat_newt, taxa_are_rows = TRUE)
TAX = tax_table(taxmat_newt)

physeq = phyloseq(OTU, TAX)


random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

physeq2 = phyloseq(OTU, TAX, sampledata[match(sampledata$X.OTU.ID,colnames(otumat_newt)),], random_tree)

unifrac_newt = UniFrac(physeq2, weighted=TRUE)
unifrac_newt = data.frame(read.csv("data/weightedunifrac.csv",row.names = 1))
unifrac_newt = unifrac_newt[as.numeric(as.character(na.exclude(match(sampledata$X.OTU.ID,rownames(unifrac_newt))))),
                            as.numeric(as.character(na.exclude(match(sampledata$X.OTU.ID,colnames(unifrac_newt)))))]

newt_NMDS=metaMDS(unifrac_newt,k=2,trymax=10000)

# from https://chrischizinski.github.io/rstats/vegan-ggplot2/
data.scores <- as.data.frame(scores(newt_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- as.character(sampledata$Disturbance[match(rownames(data.scores),sampledata$X.OTU.ID)])  #  add the grp variable created earlier

data.scores$grp <- as.character(sampledata$Season[na.exclude(match(rownames(data.scores),sampledata$X.OTU.ID))])
head(data.scores)  #look at the data

seasonpredisturb = ggplot() +  
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  scale_colour_manual(values=c("Spring" = "red", "Summer" = "blue",
                               "Fall" = "green", "Winter" = "purple"),
                      breaks = c("Spring","Summer","Fall","Winter")) +
  ylim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  xlim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  annotate("text", x = max(data.scores$NMDS1,data.scores$NMDS2) - 0.07, 
           y = min(data.scores$NMDS2,data.scores$NMDS1) + 0.01,
           label = paste0("Stress=",as.character(round(newt_NMDS$stress,3)))) + 
  scale_shape_manual(values=c("Spring" = 17, "Summer" = 15,
                              "Fall" = 16, "Winter" = 3),
                    breaks = c("Spring","Summer","Fall","Winter"))+
  theme_bw() + ylab("") + xlab("")+ 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.position = c(0.01,0.01),#c(0.97, 0.99),
        legend.justification = c("left", "bottom"),legend.background = element_rect(fill = NA))#,color="black"))
  coord_equal()


##### by season for post disturb 1##### 
otumat_newt = read.csv("./data/otutable.csv",row.names = 1)
taxmat_newt = as.matrix(read.csv("./data/taxtable.csv",row.names = 1))
sampledata = read.csv("./data/sampledata.csv")#, row.names = 1) 
sampledata = sampledata %>% dplyr::filter(Disturbance == "Post disturbance 1")
sampledata$X.OTU.ID = as.character(sampledata$X.OTU.ID)
otumat_newt = otumat_newt[,as.numeric(as.character(na.exclude(match(sampledata$X.OTU.ID,colnames(otumat_newt)))))]

OTU = otu_table(otumat_newt, taxa_are_rows = TRUE)
TAX = tax_table(taxmat_newt)

physeq = phyloseq(OTU, TAX)


random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

physeq2 = phyloseq(OTU, TAX, sampledata[match(sampledata$X.OTU.ID,colnames(otumat_newt)),], random_tree)

unifrac_newt = UniFrac(physeq2, weighted=TRUE)
unifrac_newt = data.frame(read.csv("data/weightedunifrac.csv",row.names = 1))
unifrac_newt = unifrac_newt[as.numeric(as.character(na.exclude(match(sampledata$X.OTU.ID,rownames(unifrac_newt))))),
                            as.numeric(as.character(na.exclude(match(sampledata$X.OTU.ID,colnames(unifrac_newt)))))]

newt_NMDS=metaMDS(unifrac_newt,k=2,trymax=10000)

# from https://chrischizinski.github.io/rstats/vegan-ggplot2/
data.scores <- as.data.frame(scores(newt_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- as.character(sampledata$Disturbance[match(rownames(data.scores),sampledata$X.OTU.ID)])  #  add the grp variable created earlier

data.scores$grp <- as.character(sampledata$Season[na.exclude(match(rownames(data.scores),sampledata$X.OTU.ID))])
head(data.scores)  #look at the data

seasonpostdisturb1 = ggplot() +  
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  scale_colour_manual(values=c("Spring" = "red", "Summer" = "blue",
                               "Fall" = "green", "Winter" = "purple"),
                      breaks = c("Spring","Summer","Fall","Winter")) +
  ylim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  xlim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  scale_shape_manual(values=c("Spring" = 17, "Summer" = 15,
                              "Fall" = 16, "Winter" = 3),
                     breaks = c("Spring","Summer","Fall","Winter")) +
  annotate("text", x = max(data.scores$NMDS1,data.scores$NMDS2) - 0.07, 
           y = min(data.scores$NMDS2,data.scores$NMDS1) + 0.01,
           label = paste0("Stress=",as.character(round(newt_NMDS$stress,3)))) +
  theme_bw() + ylab("") + xlab("") +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.position = c(0.01,0.01),#c(0.97, 0.99),
        legend.justification = c("left", "bottom"),legend.background = element_rect(fill = NA))

cowplot::plot_grid(disturb,seasonpredisturb,seasonpostdisturb1,ncol=1,labels = c("A","B","C"))
ggsave("nmds_ords.pdf", height = 4 * 3, width = 4.5, units = "in")

#### NMDS Jaccard ####

metrich = read.csv("data/metabolite.csv")
sampledata = read.csv("./data/sampledata.csv")
metrich["Disturbance"] = as.character(sampledata$Disturbance[match(as.character(metrich$Date), 
                                                                   as.character(sampledata$Date))])
jaccardmatx = as.matrix(dist(metrich[,1:189], method="jaccard"))
jaccardmatrix = matrix(0, nrow = nrow(metrich), ncol = nrow(metrich))
for (p in c(1:nrow(metrich))) {
  for (q in c(1:nrow(metrich))) {
  numerator = length(intersect(colnames(metrich)[which(!is.na(as.numeric(metrich[p,1:189])))], 
           colnames(metrich)[which(!is.na(as.numeric(metrich[q,1:189])))])) - 1
  
  denominator = length(union(colnames(metrich)[which(!is.na(as.numeric(metrich[p,1:189])))], 
                               colnames(metrich)[which(!is.na(as.numeric(metrich[q,1:189])))])) - 1
  jaccardmatrix[p,q] = numerator/denominator
  }
}
metrichcts = metrich[,2:189]
metrichcts[is.na(metrichcts)] <- 0
metrichcts[metrichcts > 0] <- 1
jacvegan = 1 - vegan::vegdist(metrichcts, method="jaccard")

newt_NMDS=metaMDS(metrichcts, distance = "jaccard",binary = TRUE, trace = FALSE, trymax = 1000, autotransform = FALSE)

data.scores <- as.data.frame(scores(newt_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- as.character(metrich$Disturbance)  
head(data.scores)

alldataplot = ggplot() +  
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  scale_colour_manual(name = "", values=c("Pre-disturbance" = "red", "Post disturbance 1" = "blue",
                                          "Post disturbance 2" = "green"),
                      breaks = c("Pre-disturbance","Post disturbance 1",
                                 "Post disturbance 2"),
                      labels = c("Pre-disturbance","Post disturbance 1",
                                 "Post disturbance 2")) +
  annotate("text", x = max(data.scores$NMDS1,data.scores$NMDS2) - 0.2, 
           y = min(data.scores$NMDS2,data.scores$NMDS1) + 0.01,
           label = paste0("Stress=",as.character(round(newt_NMDS$stress,3)))) +
  ylim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  xlim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  scale_shape_manual(name = "", values = c("Pre-disturbance" = 15,
                                           "Post disturbance 1" = 16,
                                           "Post disturbance 2" = 17),
                     breaks = c("Pre-disturbance","Post disturbance 1",
                                "Post disturbance 2"))+
  theme_bw() + ylab("") + xlab("")+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.position = c(0.01,0.01),#c(0.97, 0.99),
        legend.justification = c("left", "bottom"),legend.background = element_rect(fill = NA))

##### Do Jaccard matrix based on only pre-disturbance, group by season ##### 

metrich = read.csv("data/metabolite.csv")
sampledata = read.csv("./data/sampledata.csv")
metrich["Disturbance"] = as.character(sampledata$Disturbance[match(as.character(metrich$Date), as.character(sampledata$Date))])
metrich["Season"] = as.character(sampledata$Season[match(metrich$Date, sampledata$Date)])
metrich = metrich %>% dplyr::filter(Disturbance == "Pre-disturbance")
metrichcts = metrich[,2:189]
metrichcts[is.na(metrichcts)] <- 0
metrichcts[metrichcts > 0] <- 1
newt_NMDS=metaMDS(metrichcts, distance = "jaccard",binary = TRUE, trace = FALSE, trymax = 1000, autotransform = FALSE)

data.scores <- as.data.frame(scores(newt_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- as.character(metrich$Season)  
head(data.scores)

seasonpredisturb_jaccard = ggplot() +  
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  scale_colour_manual(values=c("Spring" = "red", "Summer" = "blue",
                               "Fall" = "green", "Winter" = "purple"),
                      breaks = c("Spring","Summer","Fall","Winter")) +
  scale_shape_manual(values=c("Spring" = 17, "Summer" = 15,
                              "Fall" = 16, "Winter" = 3),
                     breaks = c("Spring","Summer","Fall","Winter")) +
  annotate("text", x = max(data.scores$NMDS1,data.scores$NMDS2) - 0.2, 
           y = min(data.scores$NMDS2,data.scores$NMDS1) + 0.01,
           label = paste0("Stress=",as.character(round(newt_NMDS$stress,3)))) +
  ylim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  xlim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  theme_bw() + ylab("") + xlab("")+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.position = c(0.01,0.01),#c(0.97, 0.99),
        legend.justification = c("left", "bottom"),legend.background = element_rect(fill = NA))

##### Same thing for post disturbance 1 #####

metrich = read.csv("data/metabolite.csv")
sampledata = read.csv("./data/sampledata.csv")
metrich["Disturbance"] = as.character(sampledata$Disturbance[match(as.character(metrich$Date), as.character(sampledata$Date))])
metrich["Season"] = as.character(sampledata$Season[match(as.character(metrich$Date), as.character(sampledata$Date))])
metrich = metrich %>% dplyr::filter(Disturbance == "Post disturbance 1")
metrichcts = metrich[,2:189]
metrichcts[is.na(metrichcts)] <- 0
metrichcts[metrichcts > 0] <- 1
newt_NMDS=metaMDS(metrichcts, distance = "jaccard",binary = TRUE, trace = FALSE, trymax = 1000, autotransform = FALSE)

data.scores <- as.data.frame(scores(newt_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- as.character(metrich$Season)  
head(data.scores)

seasonpostdisturb1_jaccard = ggplot() +  
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=3) + # add the point markers
  scale_colour_manual(values=c("Spring" = "red", "Summer" = "blue",
                               "Fall" = "green", "Winter" = "purple"),
                      breaks = c("Spring","Summer","Fall","Winter")) +
  scale_shape_manual(values=c("Spring" = 17, "Summer" = 15,
                              "Fall" = 16, "Winter" = 3),
                     breaks = c("Spring","Summer","Fall","Winter")) +
  annotate("text", x = max(data.scores$NMDS1,data.scores$NMDS2) - 0.2, 
           y = min(data.scores$NMDS2,data.scores$NMDS1) + 0.01,
           label = paste0("Stress=",as.character(round(newt_NMDS$stress,3)))) +
  ylim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  xlim(c(min(data.scores$NMDS1,data.scores$NMDS2), max(data.scores$NMDS1, data.scores$NMDS2))) + 
  theme_bw() + ylab("") + xlab("")+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())+
  theme(legend.position = c(0.01,0.01),#c(0.97, 0.99),
        legend.justification = c("left", "bottom"),legend.background = element_rect(fill = NA))

cowplot::plot_grid(alldataplot,seasonpredisturb_jaccard,seasonpostdisturb1_jaccard,ncol=1,labels = c("A","B","C"))
ggsave("nmds_ords_jaccard.pdf", height = 5 * 3, width = 5, units = "in")

