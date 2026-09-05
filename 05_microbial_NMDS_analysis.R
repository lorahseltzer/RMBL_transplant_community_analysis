# Microbial NMDS analysis: RMBL transplant experiment
#
# 1. Bacteria/archaea (16S) NMDS and Figure S1
# 2. Fungal (ITS) NMDS and Figure S2

library(vegan)
library(tidyverse)
library(cowplot)
library(ggpubr)

# Bacteria/archaea NMDS --------------------------------------------------

# Load filtered 16S data.
topbact <- read.csv("data/bacteria_archaea_asv_table.csv")

# Get info I want from topbact for later joining
topbact_turf <- topbact %>%
  select(Sample, turfID, season) %>%
  unite("ID_season", turfID, season, sep=" ", remove=FALSE) %>%
  distinct(turfID, Sample, season, ID_season)

# Grouping data: Make a dataframe of how data is grouped (sample, treatmentOrigin, season) for use with plotting
bactgroup <- topbact %>%
  ungroup() %>%
  select(Sample, treatmentOrigin, season) %>% 
  distinct(Sample, .keep_all = TRUE)

# Community data: Make a matrix of sample by OTU (filled in with abundance= read counts) to use to find the distances
bactdist <- topbact %>%
  ungroup() %>% 
  select(Sample, ASV_id, Abundance) %>% 
  pivot_wider(names_from = ASV_id, values_from = Abundance,values_fill=0 ) %>% 
  column_to_rownames("Sample")
# Should I be using relative abundance instead of raw read counts? No. If you use raw read counts, the overall abundance at a site will be considered part of the difference between samples (size and shape of the count vectors will be considered). If you use the relative abundances, it will only take into account the shape differences. https://econ.upf.edu/~michael/stanford/maeb5.pdf With Bray-Curtis, 0 means the samples are the same, and 1 means they are completely different.

# NMDS ordination
# metaMDS applies a square root transformation and calculates Bray-Curtis distances for the community matrix
set.seed(1000) # set seed of random number generator used by metaMDS
bactdist2 <- metaMDS(bactdist, 
                     distance="bray",
                     k=2, # selected number of dimensions
                     maxit=999, # max iterations
                     trymax = 500,
                     wascores = TRUE) # calculates species scores)
bactdist2 # display the results
# You want stress level to be low. Higher than 0.2 is poor (risks for false interpretation). 0.1 - 0.2 is fair (some distances can be misleading for interpretation). 0.05 - 0.1 is good (can be confident in inferences from plot). Less than 0.05 is excellent (this can be rare). (from https://rpubs.com/CPEL/NMDS)

stressplot(bactdist2) # "Shepard plot, which shows scatter around the regression between the interpoint distances in the final configuration (i.e., the distances between each pair of communities) against their original dissimilarities. Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions." 
ordiplot(bactdist2, display = "sites") # basic plot of each sample
ordiplot(bactdist2, display = "species") # basic plot of each ASV

###########
## 2. Figure S1

# Extract the site scores for use in ggplot (i.e. the score for each sample)
sample.scores <-as.data.frame(scores(bactdist2, "sites"))
sample.scores$Sample <- rownames(sample.scores)

# Join metadata by sample ID rather than relying on row order.
sample.scores <- sample.scores %>%
  left_join(bactgroup, by = "Sample") %>%
  left_join(
    topbact_turf %>% select(Sample, turfID, ID_season),
    by = "Sample"
  )
sample.scores$treatmentOrigin <- factor(sample.scores$treatmentOrigin, levels = c("High Control", "Mid Control", "Low Control", "High W2", "High W1", "Mid W1", "Low C2", "Low C1", "Mid C1"))
sample.scores$season <- factor(sample.scores$season, levels = c("early", "peak", "late"))

# Extract the species scores for use in ggplot (i.e. the score for each ASV)
OTU.scores <- as.data.frame(scores(bactdist2, "species"))  #Using the scores function from vegan to extract the OTU scores and convert to a data.frame
OTU.scores$OTU <- rownames(OTU.scores)  # create a column of OTU, from the rownames of OTU.scores

#Point plot of the NMDS sample scores of: Warmed without seasons in legend including only High W2
sample.scoresHighW2 <- sample.scores %>% subset(treatmentOrigin%in%c("High Control", "High W2","Low Control"))
NMDSplotHighW2 <- ggplot(sample.scoresHighW2) +
  geom_point(data=sample.scoresHighW2, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  # geom_text(data=sample.scoresHighW2, aes(x=NMDS1, y=NMDS2, label=treatmentOrigin, color = treatmentOrigin), size=2, vjust=2, hjust=0.5)+
  scale_color_manual("Origin, Treatment", values=c("black", "gray78", "firebrick"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-0.75, 0.6) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold")) #+
#ggtitle("NMDS of 16S data (Bray-Curtis): Warmed")
NMDSplotHighW2
#ggsave("output/figures/16S_NMDS_warmed_treatmentonly_onlyHighW2.png", width=7, height=7, units="in", dpi=300 )


#Point plot of the NMDS sample scores of: Warmed without seasons in legend including only High W1
sample.scoresHW1 <- sample.scores %>% subset(treatmentOrigin%in%c("High Control", "High W1","Mid Control"))
NMDSplotHighW1 <- ggplot(sample.scoresHW1) +
  geom_point(data=sample.scoresHW1, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  #geom_text(data=sample.scoresW22, aes(x=NMDS1, y=NMDS2, label=treatmentOrigin, color = treatmentOrigin), size=2, vjust=2, hjust=0.5)+
  scale_color_manual("Origin, Treatment", values=c("black", "gray56", "darkorange"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-0.75, 0.6) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold")) #+
#ggtitle("NMDS of 16S data (Bray-Curtis): Warmed")
NMDSplotHighW1
#ggsave("output/figures/16S_NMDS_warmed_treatmentonly_onlyHighW1.png", width=7, height=7, units="in", dpi=300 )


#Point plot of the NMDS sample scores of: Warmed without seasons in legend including only Mid W1
sample.scoresMidW1 <- sample.scores %>% subset(treatmentOrigin%in%c("Mid Control", "Mid W1","Low Control"))
NMDSplotMidW1 <- ggplot(sample.scoresMidW1) +
  geom_point(data=sample.scoresMidW1, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  #geom_text(data=sample.scoresW22, aes(x=NMDS1, y=NMDS2, label=treatmentOrigin, color = treatmentOrigin), size=2, vjust=2, hjust=0.5)+
  scale_color_manual("Origin, Treatment", values=c("black", "gray78", "goldenrod1"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-0.75, 0.6) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold")) #+
#ggtitle("NMDS of 16S data (Bray-Curtis): Warmed")
NMDSplotMidW1
#ggsave("output/figures/16S_NMDS_warmed_treatmentonly_onlyMidW1.png", width=7, height=7, units="in", dpi=300 )


#Point plot of the NMDS sample scores of: Cooled
sample.scoresC <- sample.scores %>% subset(treatmentOrigin%in%c("High Control", "Low C2", "Low C1", "Low Control"))
NMDSplotC <- ggplot(sample.scoresC) +
  geom_point(data=sample.scoresC, aes(x=NMDS1, y=NMDS2, shape=season, color=treatmentOrigin), size=4) +
  #geom_text(data=sample.scoresC, aes(x=NMDS1, y=NMDS2, label=treatmentOrigin, color = treatmentOrigin), size=2, vjust=2, hjust=0.5)+
  scale_color_manual("Origin, Treatment", values=c("black", "gray78", "lightskyblue", "blue4"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold")) +
  ggtitle("NMDS of 16S data: Cooled")
NMDSplotC
#ggsave("output/figures/16S_NMDS_cooled.png", width=7, height=7, units="in", dpi=300 )


#Point plot of the NMDS sample scores of: Cooled, no season, only Low C2
sample.scoresC2 <- sample.scores %>% subset(treatmentOrigin%in%c("High Control", "Low C2", "Low Control"))
NMDSplotC2 <- ggplot(sample.scoresC2) +
  geom_point(data=sample.scoresC2, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  #geom_text(data=sample.scoresC2, aes(x=NMDS1, y=NMDS2, label=treatmentOrigin, color = treatmentOrigin), size=2, vjust=2, hjust=0.5)+
  scale_color_manual("Origin, Treatment", values=c("black", "gray78", "blue4"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-0.75,1) +
  xlim(-1.5,1.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold"))# +
#ggtitle("NMDS of 16S data (Bray-Curtis): Cooled")
NMDSplotC2
#ggsave("output/figures/16S_NMDS_cooled_onlyLowC2.png", width=7, height=7, units="in", dpi=300 )


#Point plot of the NMDS sample scores of: Cooled, no season, only Low C1
sample.scoresC1 <- sample.scores %>% subset(treatmentOrigin%in%c("Mid Control", "Low C1", "Low Control"))
NMDSplotC1 <- ggplot(sample.scoresC1) +
  geom_point(data=sample.scoresC1, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  #geom_text(data=sample.scoresC1, aes(x=NMDS1, y=NMDS2, label=treatmentOrigin, color = treatmentOrigin), size=2, vjust=2, hjust=0.5)+
  scale_color_manual("Origin, Treatment", values=c("gray56", "gray78", "dodgerblue"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-0.75,1) +
  xlim(-1.5,1.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold"))# +
#ggtitle("NMDS of 16S data (Bray-Curtis): Cooled")
NMDSplotC1
#ggsave("output/figures/16S_NMDS_cooled_onlyLowC1.png", width=7, height=7, units="in", dpi=300 )

#Point plot of the NMDS sample scores of: Cooled, no season, only Mid C1
sample.scoresMC1 <- sample.scores %>% subset(treatmentOrigin%in%c("Mid Control", "Mid C1", "High Control"))
NMDSplotMC1 <- ggplot(sample.scoresMC1) +
  geom_point(data=sample.scoresMC1, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  #geom_text(data=sample.scoresMC1, aes(x=NMDS1, y=NMDS2, label=treatmentOrigin, color = treatmentOrigin), size=2, vjust=2, hjust=0.5)+
  scale_color_manual("Origin, Treatment", values=c("black", "gray56", "lightskyblue"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-0.75,1) +
  xlim(-1.5,1.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold")) #+
#ggtitle("NMDS of 16S data (Bray-Curtis): Cooled")
NMDSplotMC1
#ggsave("output/figures/16S_NMDS_cooled_onlyMidC1.png", width=7, height=7, units="in", dpi=300 )

# Boxplots of the NMDS1 mean sample scores of: warmed
boxNMDS1W <- ggplot(sample.scores, aes(x=treatmentOrigin, y=NMDS1, color=treatmentOrigin)) +
  scale_color_manual("Origin site & treatment", values=c("black", "gray56", "gray82", "firebrick", "darkorange", "goldenrod1"))+
  geom_boxplot(alpha=0.05) +
  ylim(-0.75, 0.5) +
  xlab("Origin and Treatment") +
  scale_x_discrete(limits=c("High Control", "Mid Control", "Low Control", "High W2", "High W1", "Mid W1")) + # show only these treat
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxNMDS1W 

boxNMDS2W <- ggplot(sample.scores, aes(x=treatmentOrigin, y=NMDS2, color=treatmentOrigin)) +
  scale_color_manual("Origin site & treatment", values=c("black", "gray56", "gray82", "firebrick", "darkorange", "goldenrod1"))+
  geom_boxplot(alpha=0.05) +
  ylim(-0.75, 0.5) +
  xlab("Origin and Treatment") +
  scale_x_discrete(limits=c("High Control", "Mid Control", "Low Control", "High W2", "High W1", "Mid W1")) + # show these treat
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxNMDS2W

boxplots16SNMDS <- ggarrange(boxNMDS1W, boxNMDS2W,
                     #labels = c("A", "B"),
                     ncol = 2, nrow = 1)
boxplots16SNMDS
#ggsave("output/figures/16S_Box_NMDS_Warmed.png", figure3, width=6.5, height=4, units="in", dpi=300, device = png)

# Boxplots of the NMDS1 mean sample scores of: COOLED
boxNMDS1C <- ggplot(sample.scores, aes(x=treatmentOrigin, y=NMDS1, color=treatmentOrigin)) +
  scale_color_manual("Origin site & treatment", values=c("black", "gray56", "gray82", "blue4", "dodgerblue", "lightskyblue"))+
  geom_boxplot(alpha=0.05) +
  ylim(-0.75, 0.5) +
  xlab("Origin and Treatment") +
  scale_x_discrete(limits=c("High Control", "Mid Control", "Low Control", "Low C2", "Low C1", "Mid C1")) + # show only these treat
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxNMDS1C

boxNMDS2C <- ggplot(sample.scores, aes(x=treatmentOrigin, y=NMDS2, color=treatmentOrigin)) +
  scale_color_manual("Origin site & treatment", values=c("black", "gray56", "gray82", "blue4", "dodgerblue", "lightskyblue"))+
  geom_boxplot(alpha=0.05) +
  ylim(-0.75, 0.5)+
  xlab("Origin and Treatment") +
  scale_x_discrete(limits=c("High Control", "Mid Control", "Low Control", "Low C2", "Low C1", "Mid C1")) + # show these treat
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxNMDS2C

figure4 <- ggarrange(boxNMDS1C, boxNMDS2C,
                     #labels = c("A", "B"),
                     ncol = 2, nrow = 1)
figure4

# Assemble Figure S1 after all component panels exist.
BacteriaNMDS_Warm <- plot_grid(
  NMDSplotHighW2, NMDSplotHighW1, NMDSplotMidW1, boxplots16SNMDS,
  labels = c("A", "B", "C", "D"), label_size = 20
)
BacteriaNMDS_Cool <- plot_grid(
  NMDSplotC2, NMDSplotC1, NMDSplotMC1, figure4,
  labels = c("E", "F", "G", "H"), label_size = 20
)
figure_s1 <- plot_grid(BacteriaNMDS_Warm, BacteriaNMDS_Cool, ncol = 1)
ggsave(
  "output/figures/figure_s1_bacteria_archaea_nmds.png",
  plot = figure_s1, width = 12, height = 16, units = "in", dpi = 300
)




# Fungal NMDS ------------------------------------------------------------

topfungi <- read.csv("data/fungal_asv_table.csv")
# Get info I want from topfungi for later joining
topfungi_turf <- topfungi %>%
  select(Sample, turfID, season) %>%
  unite("ID_season", turfID, season, sep=" ", remove=FALSE) %>%
  distinct(turfID, Sample, season, ID_season)

# Grouping data: Make a dataframe of how data is grouped (sample, treatmentOrigin, season) for use with plotting
fungigroup <- topfungi %>%
  ungroup() %>%
  select(Sample, treatmentOrigin, season) %>% 
  distinct(Sample, .keep_all = TRUE)

# Community data: Make a matrix of sample by OTU (filled in with abundance = read counts) to use to find the distances
fungidist <- topfungi %>%
  ungroup()%>%
  select(Sample, OTU, Abundance) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance,values_fill=0 ) %>% 
  column_to_rownames("Sample")

# Run NMDS ordination
set.seed(1575) # set seed of random number generator used by metaMDS
fungidist2 <- metaMDS(fungidist, 
                      distance="bray",
                      k=2, # selected number of dimensions: 2 since making 2D figures
                      maxit=999, # max iterations
                      trymax = 500,
                      wascores = TRUE) # calculates species scores
fungidist2 # display the results
# You want stress level to be low. Higher than 0.2 is poor (risks for false interpretation). 0.1 - 0.2 is fair (some distances can be misleading for interpretation). 0.05 - 0.1 is good (can be confident in inferences from plot). Less than 0.05 is excellent (this can be rare). (from https://rpubs.com/CPEL/NMDS))
stressplot(fungidist2) # "Shepard plot, which shows scatter around the regression between the interpoint distances in the final configuration (i.e., the distances between each pair of communities) against their original dissimilarities. Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions." 
ordiplot(fungidist2, display = "sites") # basic plot of each sample
ordiplot(fungidist2, display = "species") # basic plot of each OTU

# Extract the site scores for use in ggplot (i.e. the score for each sample)
sample.scoresfungi <-as.data.frame(scores(fungidist2, "sites"))
#sample.scoresfungi <-as.data.frame(fungidist2$points)
#sample.scoresfungi <- rename(sample.scoresfungi, NMDS1 = MDS1)
#sample.scoresfungi <- rename(sample.scoresfungi, NMDS2 = MDS2)
sample.scoresfungi$Sample <- rownames(sample.scoresfungi) # Add a column with sample names.
# Join metadata by sample ID rather than relying on row order.
sample.scoresfungi <- sample.scoresfungi %>%
  left_join(fungigroup, by = "Sample") %>%
  left_join(
    topfungi_turf %>% select(Sample, turfID, ID_season),
    by = "Sample"
  )
sample.scoresfungi$treatmentOrigin <- factor(sample.scoresfungi$treatmentOrigin, levels = c("High Control", "Mid Control", "Low Control", "High W2", "High W1", "Mid W1", "Low C2", "Low C1", "Mid C1"))
sample.scoresfungi$season <- factor(sample.scoresfungi$season, levels = c("early", "peak", "late"))

# Extract the species scores for use in ggplot (i.e. the score for each ASV)
OTU.scoresfungi <- as.data.frame(scores(fungidist2, "species"))  #Using the scores function from vegan to extract the OTU scores and convert to a data.frame
OTU.scoresfungi$OTU <- rownames(OTU.scoresfungi)  # create a column of OTU, from the rownames of OTU.scores




#### 4. Figure S2
#NMDS Warmed (without season, only W2)
sample.scoresWf2 <- sample.scoresfungi %>% subset(treatmentOrigin%in%c("High Control", "Low Control", "High W2"))
NMDSplotWf2 <- ggplot(sample.scoresWf2) +
  geom_point(data=sample.scoresWf2, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  scale_color_manual("Origin, Treatment", values=c("black", "gray78", "firebrick"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-1.25, 1.25) +
  xlim(-1.25,1.25)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold")) 
NMDSplotWf2

#NMDS Warmed (without season, only High W1)
sample.scoresWf3 <- sample.scoresfungi %>% subset(treatmentOrigin%in%c("High Control", "Mid Control", "High W1"))
NMDSplotWf3 <- ggplot(sample.scoresWf3) +
  geom_point(data=sample.scoresWf3, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  scale_color_manual("Origin, Treatment", values=c("black", "gray56", "darkorange"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-1.25, 1.25) +
  xlim(-1.25,1.25)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold")) 
NMDSplotWf3

#NMDS Warmed (without season, only Mid W1)
sample.scoresWf4 <- sample.scoresfungi %>% subset(treatmentOrigin%in%c("Mid Control", "Low Control", "Mid W1"))
NMDSplotWf4 <- ggplot(sample.scoresWf4) +
  geom_point(data=sample.scoresWf4, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  scale_color_manual("Origin, Treatment", values=c("gray56", "gray78", "goldenrod1"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-1.25, 1.25) +
  xlim(-1.25,1.25) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold"))
NMDSplotWf4




#NMDS Cooled (with season in legend, treatments from low)
sample.scoresCf <- sample.scoresfungi %>% subset(treatmentOrigin%in%c("High Control", "Low C2", "Low C1", "Low Control"))
NMDSplotCf <- ggplot(sample.scoresCf) +
  geom_point(data=sample.scoresCf, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  scale_color_manual("Origin, Treatment", values=c("black", "gray78", "lightskyblue", "blue4"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ggtitle("NMDS of ITS data (Bray-Curtis): Cooled") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold")) 
NMDSplotCf

# NMDS Cooled (no season, only Low C2)
sample.scoresCf2 <- sample.scoresfungi %>% subset(treatmentOrigin%in%c("High Control", "Low Control", "Low C2"))
NMDSplotCf2 <- ggplot(sample.scoresCf2) +
  geom_point(data=sample.scoresCf2, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  scale_color_manual("Origin, Treatment", values=c("black", "gray78", "blue4"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-1.25,1) +
  xlim(-1.25,1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold"))
NMDSplotCf2

# NMDS Cooled (no season, only Low C1)
sample.scoresCf3 <- sample.scoresfungi %>% subset(treatmentOrigin%in%c("Low Control", "Low C1", "Mid Control"))
NMDSplotCf3 <- ggplot(sample.scoresCf3) +
  geom_point(data=sample.scoresCf3, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  scale_color_manual("Origin, Treatment", values=c("gray56", "gray78", "dodgerblue"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-1.25,1) +
  xlim(-1.25,1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold"))
NMDSplotCf3


# NMDS Cooled (no season, only Mid C1)
sample.scoresCf4 <- sample.scoresfungi %>% subset(treatmentOrigin%in%c("Mid Control", "Mid C1", "High Control"))
NMDSplotCf4 <- ggplot(sample.scoresCf4) +
  geom_point(data=sample.scoresCf4, aes(x=NMDS1, y=NMDS2, color=treatmentOrigin), size=4) +
  scale_color_manual("Origin, Treatment", values=c("black", "gray56", "lightskyblue"))+
  stat_ellipse(aes(x=NMDS1, y=NMDS2, color=treatmentOrigin)) +
  ylim(-1.25,1) +
  xlim(-1.25,1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.title= element_text(size = 15), 
        axis.text = element_text(size=10), 
        strip.text=element_text(size=15), 
        plot.title = element_text(size=22), 
        legend.text=element_text(size=15), 
        legend.title=element_text(size=15, face="bold"))
NMDSplotCf4


# Boxplots of the NMDS1 mean sample scores of: warmed
boxNMDS1Wf <- ggplot(sample.scoresfungi, aes(x=treatmentOrigin, y=NMDS1, color=treatmentOrigin)) +
  scale_color_manual("Origin site & treatment", values=c("black", "gray56", "gray82", "firebrick", "darkorange", "goldenrod1"))+
  geom_boxplot(alpha=0.05) +
  ylim(-0.6, 0.8)+
  xlab("Origin and Treatment") +
  scale_x_discrete(limits=c("High Control", "Mid Control", "Low Control", "High W2", "High W1", "Mid W1")) + # show only these treat
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxNMDS1Wf 

boxNMDS2Wf <- ggplot(sample.scoresfungi, aes(x=treatmentOrigin, y=NMDS2, color=treatmentOrigin)) +
  scale_color_manual("Origin site & treatment", values=c("black", "gray56", "gray82", "firebrick", "darkorange", "goldenrod1"))+
  geom_boxplot(alpha=0.05) +
  ylim(-0.6, 0.8)+
  xlab("Origin and Treatment") +
  scale_x_discrete(limits=c("High Control", "Mid Control", "Low Control", "High W2", "High W1", "Mid W1")) + # show these treat
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxNMDS2Wf

figure3f <- ggarrange(boxNMDS1Wf, boxNMDS2Wf,
                      #labels = c("A", "B"),
                      ncol = 2, nrow = 1)
figure3f


# Boxplots of the NMDS1 mean sample scores of: COOLED
boxNMDS1Cf <- ggplot(sample.scoresfungi, aes(x=treatmentOrigin, y=NMDS1, color=treatmentOrigin)) +
  scale_color_manual("Origin site & treatment", values=c("black", "gray56", "gray82", "blue4", "dodgerblue", "lightskyblue"))+
  geom_boxplot(alpha=0.05) +
  xlab("Origin and Treatment") +
  scale_x_discrete(limits=c("High Control", "Mid Control", "Low Control", "Low C2", "Low C1", "Mid C1")) + # show only these treat
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxNMDS1Cf

boxNMDS2Cf <- ggplot(sample.scoresfungi, aes(x=treatmentOrigin, y=NMDS2, color=treatmentOrigin)) +
  scale_color_manual("Origin site & treatment", values=c("black", "gray56", "gray82", "blue4", "dodgerblue", "lightskyblue"))+
  geom_boxplot(alpha=0.05) +
  xlab("Origin and Treatment") +
  scale_x_discrete(limits=c("High Control", "Mid Control", "Low Control", "Low C2", "Low C1", "Mid C1")) + # show these treat
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
boxNMDS2Cf

figure4f <- ggarrange(boxNMDS1Cf, boxNMDS2Cf,
                      #labels = c("A", "B"),
                      ncol = 2, nrow = 1)
figure4f

# Assemble Figure S2 after all component panels exist.
FungiNMDS_Warm <- plot_grid(
  NMDSplotWf2, NMDSplotWf3, NMDSplotWf4, figure3f,
  labels = c("A", "B", "C", "D"), label_size = 20
)
FungiNMDS_Cool <- plot_grid(
  NMDSplotCf2, NMDSplotCf3, NMDSplotCf4, figure4f,
  labels = c("E", "F", "G", "H"), label_size = 20
)
figure_s2 <- plot_grid(FungiNMDS_Warm, FungiNMDS_Cool, ncol = 1)
ggsave(
  "output/figures/figure_s2_fungi_nmds.png",
  plot = figure_s2, width = 12, height = 16, units = "in", dpi = 300
)

# Record NMDS stress values used to assess ordination fit.
nmds_stress <- data.frame(
  community = c("bacteria_archaea", "fungi"),
  stress = c(bactdist2$stress, fungidist2$stress)
)
write.csv(
  nmds_stress,
  "output/tables/microbial_nmds_stress.csv",
  row.names = FALSE
)
