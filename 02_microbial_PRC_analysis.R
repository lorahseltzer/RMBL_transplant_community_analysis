## 02_microbial_PRC_analysis.R
# Microbial community analysis: RMBL transplant experiment
#
# 1. Bacteria/archaea (16S) Principal Response Curve analysis
# 2. Fungal (ITS) Principal Response Curve analysis
# 3. Species-weight correlations and Figure 3


library(vegan)
library(ggvegan)
library(ggplot2)
library(tidyverse)
library(cowplot)

# Make permutation tests reproducible.
set.seed(20250902)

# Bacteria/archaea data preparation --------------------------------------

# Load filtered 16S data.
topbact <- read.csv("data/bacteria_archaea_asv_table.csv")

# Correct site elevation from 3200 to 3150 m.
topbact <- topbact %>% mutate(destElevm = ifelse(destElevm == 3200, 3150, destElevm)) %>% mutate(originElevm = ifelse(originElevm == 3200, 3150, originElevm))

# Retain sample metadata for later joins.
topbactInfo <- topbact %>%
  select(Sample, turfID, season, destElevm, originElevm, treatmentOrigin) %>%
  unite("ID_season", turfID, season, sep=" ", remove=FALSE) %>%
  distinct(turfID, Sample, season, ID_season, destElevm, originElevm, treatmentOrigin)

# Log-transform read abundance and encode seasonal order.
topbactASV <- topbact %>%
  mutate(logAbundance= log(Abundance)) %>%
  filter(!is.na(ASV_id), !is.na(logAbundance), !is.na(season), !is.na(turfID)) %>%
  mutate(seasonnum = season) %>%
  mutate(seasonnum= recode(seasonnum, "late"="1", "early"="2", "peak"="3")) %>%
  select(turfID, season, treatmentOrigin, ASV_id, logAbundance, seasonnum)


# Bacteria/archaea PRCs ---------------------------------------------------
# Warming treatments.
bactFat_Warming <- topbactASV %>% 
  filter(treatmentOrigin %in% c("High Control", "High W2", "High W1", "Mid W1", "Low Control", "Mid Control")) %>%
  mutate(ASV_id = paste0("ASV_", ASV_id)) %>%
  spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
  mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels = c("High Control", "High W2", "High W1","Mid W1", "Mid Control", "Low Control"))) %>%
  ungroup() 

# Retain only numeric community-response columns.
bactData_Warming <- bactFat_Warming %>%
  ungroup() %>%
  select(-turfID, -season, -treatmentOrigin, -seasonnum)

fit_bactWarming <- prc(response = bactData_Warming, treatment = bactFat_Warming$treatmentOrigin, time = bactFat_Warming$seasonnum, scale=TRUE)

summary(fit_bactWarming)
plot(fit_bactWarming, scaling="species", ylab="Warming effect on proportional \n dominance of plant species", xlab= "seasonal time point")


# Cooling treatments.
bactFat_Cooling <- topbactASV %>% 
  filter(treatmentOrigin %in% c("Low Control", "Low C2", "Low C1", "Mid C1", "Mid Control", "High Control")) %>%
  mutate(ASV_id = paste0("ASV_", ASV_id)) %>%
  spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
  mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels = c("Low Control", "Low C2", "Low C1", "Mid C1", "Mid Control", "High Control"))) %>%
  ungroup()

# Retain only numeric community-response columns.
bactData_Cooling <- bactFat_Cooling %>%
  ungroup() %>%
  select(-turfID, -season, -treatmentOrigin, -seasonnum)

fit_bactCooling <- prc(response = bactData_Cooling, treatment = bactFat_Cooling$treatmentOrigin, time = bactFat_Cooling$seasonnum, scale=TRUE)

summary(fit_bactCooling)
plot(fit_bactCooling, scaling="species", ylab="Cooling effect on proportional \n dominance of plant species", xlab= "Year")

# Bacteria/archaea species weights ---------------------------------------

# Extract warming species weights.
bactscoresdf <- fortify(fit_bactWarming) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="warmed_bacteria_gradient") %>%
  rownames_to_column(var="species") %>%
  mutate(species = as.numeric(sub("^ASV_", "", species)))

# Extract cooling species weights.
bactscoresdfC <- fortify(fit_bactCooling) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="cooled_bacteria_gradient") %>%
  rownames_to_column(var="species") %>%
  mutate(species = as.numeric(sub("^ASV_", "", species)))

# Combine
bactloading <- bind_rows(bactscoresdf, bactscoresdfC) %>% 
  select(species, Response, community) %>%
  dplyr::rename("ASV_id" = "species") %>%
  pivot_wider(id_cols=ASV_id, names_from=community, values_from=Response)
bactloading$ASV_id <- as.numeric(bactloading$ASV_id)

topbactTaxinfo <- topbact %>% distinct(ASV_id, .keep_all=TRUE) %>% select(ASV_id, Genus, Family, Class, Order, Phylum)
bactloading <- left_join(bactloading, topbactTaxinfo, by="ASV_id")



# Bacteria/archaea species-weight correlations ---------------------------

# Compare experimental warming and cooling weights with weights along the
# elevation gradient, and test for symmetry between warming and cooling.

bactW_origin <- topbactASV %>% 
  filter(treatmentOrigin %in% c("High Control", "High W2", "High W1", "Mid W1")) %>%
  mutate(ASV_id = paste0("ASV_", ASV_id)) %>%
  spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
  mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels = c("High Control", "High W2", "High W1","Mid W1"))) %>%
  ungroup() 
# Make object containing only numerical cover data. 
bactWdata_origin<- bactW_origin %>%
  ungroup() %>%
  select(-turfID, -season, -treatmentOrigin, -seasonnum)
fit_bactW_origin <- prc(response = bactWdata_origin, treatment = bactW_origin$treatmentOrigin, time = bactW_origin$seasonnum, scale=TRUE)
# get bact species weights 
bactscoresW_or<- fortify(fit_bactW_origin) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="warmed_bacteria") %>%
  rownames_to_column(var="species") %>%
  mutate(species = as.numeric(sub("^ASV_", "", species)))


bactC_origin <- topbactASV %>% 
  filter(treatmentOrigin %in% c("Low Control", "Low C2", "Low C1", "Mid C1")) %>%
  mutate(ASV_id = paste0("ASV_", ASV_id)) %>%
  spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
  mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels = c("Low Control", "Low C2", "Low C1", "Mid C1"))) %>%
  ungroup() 
# Make object containing only numerical cover data. 
bactCdata_origin<- bactC_origin %>%
  ungroup() %>%
  select(-turfID, -season, -treatmentOrigin, -seasonnum)
fit_bactC_origin <- prc(response = bactCdata_origin, treatment = bactC_origin$treatmentOrigin, time = bactC_origin$seasonnum, scale=TRUE)
# get bact species weights 
bactscoresC_or<- fortify(fit_bactC_origin) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="cooled_bacteria") %>%
  rownames_to_column(var="species") %>%
  mutate(species = as.numeric(sub("^ASV_", "", species)))

# Test for symmetry between warming and cooling responses.
weight_wc<- bind_rows(bactscoresC_or, bactscoresW_or)
weight_wc <- weight_wc %>%
  select(species, Response, community) %>%
  dplyr::rename(coverResponse = Response) %>% 
  pivot_wider(names_from = community, values_from = coverResponse)
cor.test(weight_wc$cooled_bacteria, weight_wc$warmed_bacteria)
plot(weight_wc$cooled_bacteria, weight_wc$warmed_bacteria)

# Compare experimental warming with the natural warming gradient.
weight_wgrad <- bind_rows(bactscoresdf, bactscoresW_or)
weight_wgrad  <- weight_wgrad %>%
  select(species, Response, community) %>%
  dplyr::rename(coverResponse = Response) %>%
  pivot_wider(names_from = community, values_from = coverResponse)
cor.test(weight_wgrad $warmed_bacteria, weight_wgrad $warmed_bacteria_gradient)
plot(weight_wgrad $warmed_bacteria, weight_wgrad $warmed_bacteria_gradient)

# Compare experimental cooling with the natural cooling gradient.
weight_cgrad  <- bind_rows(bactscoresdfC, bactscoresC_or)
weight_cgrad$Response <- weight_cgrad$Response * (-1) # to reverse signs of weights to indicate flipped axis
weight_cgrad <- weight_cgrad %>%
  select(species, Response, community) %>%
  dplyr::rename(coverResponse = Response) %>%
  pivot_wider(names_from = community, values_from = coverResponse)
cor.test(weight_cgrad$cooled_bacteria, weight_cgrad$cooled_bacteria_gradient)
plot(weight_cgrad$cooled_bacteria, weight_cgrad$cooled_bacteria_gradient)

weight_bactall <- full_join(weight_cgrad, weight_wgrad, by="species", relationship="many-to-many") %>%
  select(species, warmed_bacteria_gradient, warmed_bacteria, cooled_bacteria_gradient, cooled_bacteria)
weight_bactall  <- weight_bactall  %>% mutate(species = factor(species))
write.csv(weight_bactall, "output/tables/bacteria_archaea_prc_species_weight_correlations.csv", row.names = FALSE)



# Bacteria/archaea descriptive summaries --------------------------------

# Calculate mean log-transformed abundance for each ASV in each control.
meandombact <- topbactASV %>% 
  filter(treatmentOrigin %in% c("High Control", "Low Control", "Mid Control")) %>%
  filter(!is.na(logAbundance)) %>%
  group_by(ASV_id, treatmentOrigin) %>%
  summarize(mean_abund=mean(logAbundance)) %>%
  ungroup()%>%
  mutate_if(is.numeric, round, 2) %>%
  pivot_wider(id_cols = ASV_id, names_from=treatmentOrigin, values_from=mean_abund)
meandombact <- meandombact %>% mutate(ASV_id = factor(ASV_id))
write.csv(meandombact, file = "output/tables/bacteria_archaea_mean_log_abundance_by_control.csv", row.names = FALSE)

# Join taxonomic information to the species weights.
topbactfam <- topbact %>% select(ASV_id, Genus, Family, Phylum, Class, Order) %>% distinct(ASV_id, .keep_all=TRUE)
topbactfam <- topbactfam %>% mutate(ASV_id = factor(ASV_id))

bactloading2 <- left_join(weight_bactall, meandombact, by=c("species" = "ASV_id")) %>% dplyr::rename("ASV_id"="species")
bactloading3 <- left_join(bactloading2, topbactfam, by="ASV_id") %>% mutate_if(is.numeric, round, 2)
write.csv(bactloading3, file="output/tables/bacteria_archaea_prc_species_scores_with_taxonomy.csv", row.names = FALSE)




# Fungal data preparation ------------------------------------------------

# Load filtered ITS data.
topfungi <- read.csv("data/fungal_asv_table.csv")

# Correct site elevation from 3200 to 3150 m.
topfungi <- topfungi %>% mutate(destElevm = ifelse(destElevm == 3200, 3150, destElevm)) %>% mutate(originElevm = ifelse(originElevm == 3200, 3150, originElevm))

# Retain sample metadata for later joins.
topfungiInfo <- topfungi %>%
  select(Sample, turfID, season, destElevm, originElevm, treatmentOrigin) %>%
  unite("ID_season", turfID, season, sep=" ", remove=FALSE) %>%
  distinct(turfID, Sample, season, ID_season, destElevm, originElevm, treatmentOrigin)

# Log-transform read abundance and encode seasonal order.
topfungiASV <- topfungi %>%
  mutate(logAbundance= log(Abundance)) %>%
  filter(!is.na(ASV_id), !is.na(logAbundance), !is.na(season), !is.na(turfID)) %>%
  mutate(seasonnum = season) %>%
  mutate(seasonnum= recode(seasonnum, "late"="1", "early"="2", "peak"="3")) %>%
  select(turfID, season, treatmentOrigin, ASV_id, logAbundance, seasonnum)

# Fungal PRCs ------------------------------------------------------------
# Warming treatments.
fungiFat_Warming <- topfungiASV %>% 
  filter(treatmentOrigin %in% c("High Control", "High W2", "High W1", "Mid W1", "Low Control", "Mid Control")) %>%
  mutate(ASV_id = paste0("ASV_", ASV_id)) %>%
  spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
  mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels = c("High Control", "High W2", "High W1","Mid W1", "Mid Control", "Low Control"))) %>%
  ungroup()

# Retain only numeric community-response columns.
fungiData_Warming <- fungiFat_Warming %>%
  ungroup() %>%
  select(-turfID, -season, -treatmentOrigin, -seasonnum)

fit_fungiWarming <- prc(response = fungiData_Warming, treatment = fungiFat_Warming$treatmentOrigin, time = fungiFat_Warming$seasonnum, scale=TRUE)

summary(fit_fungiWarming)
plot(fit_fungiWarming, scaling="species", ylab="Warming effect on the abundance of fungal species", xlab= "seasonal time point")

# Cooling treatments.
fungiFat_Cooling <- topfungiASV %>% 
  filter(treatmentOrigin %in% c("Low Control", "Low C2", "Low C1", "Mid C1", "Mid Control", "High Control")) %>%
  mutate(ASV_id = paste0("ASV_", ASV_id)) %>%
  spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
  mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels = c("Low Control", "Low C2", "Low C1", "Mid C1", "Mid Control", "High Control"))) %>%
  ungroup()

# Retain only numeric community-response columns.
fungiData_Cooling <- fungiFat_Cooling %>%
  ungroup() %>%
  select(-turfID, -season, -treatmentOrigin, -seasonnum)

fit_fungiCooling <- prc(response = fungiData_Cooling, treatment = fungiFat_Cooling$treatmentOrigin, time = fungiFat_Cooling$seasonnum, scale=TRUE)

summary(fit_fungiCooling)
plot(fit_fungiCooling, scaling="species", ylab="Cooling effect on abundance of fungi", xlab= "seasonal time point")


# Fungal species weights -------------------------------------------------

# Extract warming species weights.
fungiscoresdf <- fortify(fit_fungiWarming) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="warmed_fungi_gradient") %>%
  rownames_to_column(var="species") %>%
  mutate(species = as.numeric(sub("^ASV_", "", species)))

# Extract cooling species weights.
fungiscoresdfC <- fortify(fit_fungiCooling) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="cooled_fungi_gradient") %>%
  rownames_to_column(var="species") %>%
  mutate(species = as.numeric(sub("^ASV_", "", species)))

# Combine
fungiloading <- bind_rows(fungiscoresdf, fungiscoresdfC) %>% 
  select(species, Response, community) %>%
  dplyr::rename("ASV_id" = "species") %>%
  pivot_wider(id_cols=ASV_id, names_from=community, values_from=Response)
fungiloading$ASV_id <- as.numeric(fungiloading$ASV_id)

topfungiTaxinfo <- topfungi %>%
  distinct(ASV_id, .keep_all = TRUE) %>%
  select(ASV_id, Genus, Family, Class, Order, Phylum)
fungiloading <- left_join(fungiloading, topfungiTaxinfo, by="ASV_id")


# Fungal species-weight correlations ------------------------------------

# Compare experimental warming and cooling weights with weights along the
# elevation gradient, and test for symmetry between warming and cooling.

funW_origin <- topfungiASV %>% 
  filter(treatmentOrigin %in% c("High Control", "High W2", "High W1", "Mid W1")) %>%
  mutate(ASV_id = paste0("ASV_", ASV_id)) %>%
  spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
  mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels = c("High Control", "High W2", "High W1","Mid W1"))) %>%
  ungroup() 
# Make object containing only numerical cover data. 
funWdata_origin<- funW_origin %>%
  ungroup() %>%
  select(-turfID, -season, -treatmentOrigin, -seasonnum)
fit_funW_origin <- prc(response = funWdata_origin, treatment = funW_origin$treatmentOrigin, time = funW_origin$seasonnum, scale=TRUE)
# get species weights 
funscoresW_or<- fortify(fit_funW_origin) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="warmed_fungi") %>%
  rownames_to_column(var="species") %>%
  mutate(species = as.numeric(sub("^ASV_", "", species)))

funC_origin <- topfungiASV %>% 
  filter(treatmentOrigin %in% c("Low Control", "Low C2", "Low C1", "Mid C1")) %>%
  mutate(ASV_id = paste0("ASV_", ASV_id)) %>%
  spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
  mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels = c("Low Control", "Low C2", "Low C1", "Mid C1"))) %>%
  ungroup() 
# Make object containing only numerical cover data. 
funCdata_origin<- funC_origin %>%
  ungroup() %>%
  select(-turfID, -season, -treatmentOrigin, -seasonnum)
fit_funC_origin <- prc(response = funCdata_origin, treatment = funC_origin$treatmentOrigin, time = funC_origin$seasonnum, scale=TRUE)
# Extract cooling species weights.
funscoresC_or<- fortify(fit_funC_origin) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="cooled_fungi") %>%
  rownames_to_column(var="species") %>%
  mutate(species = as.numeric(sub("^ASV_", "", species)))


# Test for symmetry between warming and cooling responses.
fweight_wc<- bind_rows(funscoresC_or, funscoresW_or)
fweight_wc <- fweight_wc %>%
  select(species, Response, community) %>%
  dplyr::rename(coverResponse = Response) %>% 
  pivot_wider(names_from = community, values_from = coverResponse)
cor.test(fweight_wc$cooled_fungi, fweight_wc$warmed_fungi)
plot(fweight_wc$cooled_fungi, fweight_wc$warmed_fungi)

# Compare experimental warming with the natural warming gradient.
fweight_wgrad <- bind_rows(fungiscoresdf, funscoresW_or)
fweight_wgrad  <- fweight_wgrad %>%
  select(species, Response, community) %>%
  dplyr::rename(coverResponse = Response) %>%
  pivot_wider(names_from = community, values_from = coverResponse)
cor.test(fweight_wgrad$warmed_fungi, fweight_wgrad$warmed_fungi_gradient)
plot(fweight_wgrad$warmed_fungi, fweight_wgrad$warmed_fungi_gradient)

# Compare experimental cooling with the natural cooling gradient.
fweight_cgrad  <- bind_rows(fungiscoresdfC, funscoresC_or)
fweight_cgrad$Response <- fweight_cgrad$Response * (-1) # to reverse signs of weights to indicate flipped axis
fweight_cgrad <- fweight_cgrad %>%
  select(species, Response, community) %>%
  dplyr::rename(coverResponse = Response) %>%
  pivot_wider(names_from = community, values_from = coverResponse)
cor.test(fweight_cgrad$cooled_fungi, fweight_cgrad$cooled_fungi_gradient)
plot(fweight_cgrad$cooled_fungi, fweight_cgrad$cooled_fungi_gradient)

fweight_fungiall <- full_join(fweight_cgrad, fweight_wgrad, by="species", relationship="many-to-many") %>%
  select(species, warmed_fungi_gradient, warmed_fungi, cooled_fungi_gradient, cooled_fungi)
fweight_fungiall <- fweight_fungiall %>%
  dplyr::rename("ASV_id" = "species") %>%
  mutate(ASV_id = as.numeric(ASV_id))



# Fungal descriptive summaries ------------------------------------------
# Calculate mean log-transformed abundance for each ASV in each control.
meandomfungi <- topfungiASV %>% 
  filter(treatmentOrigin %in% c("High Control", "Low Control", "Mid Control")) %>%
  filter(!is.na(logAbundance)) %>%
  group_by(ASV_id, treatmentOrigin) %>%
  summarize(mean_abund=mean(logAbundance)) %>%
  ungroup()%>%
  mutate_if(is.numeric, round, 2) %>%
  pivot_wider(id_cols = ASV_id, names_from=treatmentOrigin, values_from=mean_abund)
meandomfungi <- meandomfungi %>% mutate(ASV_id = as.numeric(ASV_id))
write.csv(meandomfungi, file = "output/tables/fungi_mean_log_abundance_by_control.csv", row.names = FALSE)

fungiloading2 <- left_join(fweight_fungiall, meandomfungi, by="ASV_id") %>%
  left_join(topfungiTaxinfo, by = "ASV_id")
write.csv(fungiloading2, file="output/tables/fungi_prc_species_weight_correlations.csv", row.names = FALSE)


# Figure 3 plotting function --------------------------------------------
autoplot.prcWithoutSPBact <- function(object, select, xlab, ylab,
                                      title = NULL, subtitle = NULL, caption = NULL,
                                      legend.position = "top", ...) {
  ## fortify the model object
  fobj <- fortify(object, ...)
  
  ## levels of factors - do this now before we convert things
  TimeLevs <- levels(fobj$Time)
  TreatLevs <- levels(fobj$Treatment)
  
  ## convert Time to a numeric
  fobj$Time <- as.numeric(as.character(fobj$Time))
  
  ## process select
  ind <- fobj$score != "Sample"
  if(missing(select)) {
    select <- rep(TRUE,sum(ind))
  } else {
    stopifnot(isTRUE(all.equal(length(select), sum(ind))))
  }
  
  ## samples and species "scores"
  samp <- fobj[!ind, ] 
  spp <- fobj[ind,][select, ]
  
  ## base plot
  plt <- ggplot(data = samp,
                aes_string(x = 'Time', y = 'Response', group = 'Treatment',
                           colour = 'Treatment', linetype = 'Treatment'))
  ## add the control
  plt <- plt + geom_hline(yintercept = 0)
  
  ## add the coefficients
  plt <- plt + geom_line(size = 1.5) +
    theme(legend.position = legend.position, 
          legend.title = element_blank(),
          text = element_text(size=20),
          axis.text = element_text(size = 20)) +
    scale_x_continuous(breaks = c(1,2,3), labels=c("1"="2018\nlate", "2"="2019\nearly", "3"="2019\npeak"), minor_breaks = NULL)
  
  ## add labels
  if(missing(xlab)) {
    xlab <- 'Time'
  }
  if(missing(ylab)) {
    ylab <- 'Treatment'
  }
  plt <- plt + labs(x = xlab, y = ylab, title = title, subtitle = subtitle,
                    caption = caption)
  
  ## return
  plt
}

# Figure 3 ---------------------------------------------------------------
# Panel A: bacteria/archaea warming.
p2swBact <- autoplot.prcWithoutSPBact(fit_bactWarming, xlab = "Growing season time point", ylab = "Warming effect on the \nabundance of bacteria and archaea") +
  scale_colour_manual(values = c("red", "orange", "red", "orange", "red")) +
  scale_linetype_manual(values = c("dashed", "dashed", "dotted", "solid", "solid")) +
  #scale_y_continuous(limits = c(-0.3, 0.35), breaks = seq(-0.2,0.2, by=0.2)) +
  theme(
    legend.position = c(0.2, 0.8), 
    legend.title = element_blank(), 
    panel.background=element_blank(), 
    axis.line = element_line(colour = "black"), 
    legend.key.width=unit(4,"line"), 
    legend.key = element_rect(colour = "transparent"),
    axis.text=element_text(size=15))
p2swBact

# Panel B: bacteria/archaea cooling.
p2scBact <- autoplot.prcWithoutSPBact(fit_bactCooling, xlab = "Growing season time point", ylab = "Cooling effect on the \nabundance of bacteria and archaea") +
  scale_colour_manual(values = c("navyblue", "lightblue", "navy", "lightblue", "navy")) +
  scale_linetype_manual(values = c("dashed", "dashed", "dotted", "solid", "solid")) +
  #scale_y_continuous(limits = c(-0.3, 0.35), breaks = seq(-0.2,0.2, by=0.2)) +
  theme(
    legend.position = c(0.2, 0.12), 
    legend.title = element_blank(), 
    panel.background=element_blank(), 
    axis.line = element_line(colour = "black"), 
    legend.key.width=unit(4,"line"), 
    legend.key = element_rect(colour = "transparent"),
    axis.text=element_text(size=15))
p2scBact

# Combine bacteria/archaea panels A and B for interactive inspection.
RDA_speciesbact <- plot_grid(p2swBact, p2scBact, align= c("hv"), labels=c("A", "B"), label_size = 20)
RDA_speciesbact


# Panel C: fungal warming.
p2swFungi <- autoplot.prcWithoutSPBact(fit_fungiWarming, xlab = "Growing season time point", ylab = "Warming effect on the\nabundance of fungi") +
  scale_colour_manual(values = c("red", "orange", "red", "orange", "red")) +
  scale_linetype_manual(values = c("dashed", "dashed", "dotted", "solid", "solid")) +
  #scale_y_continuous(limits = c(-0.2, 0.5), breaks = seq(0,0.5, by=0.2)) +
  theme(
    legend.position = c(0.2, 0.8), 
    legend.title = element_blank(), 
    panel.background=element_blank(), 
    axis.line = element_line(colour = "black"), 
    legend.key.width=unit(4,"line"), 
    legend.key = element_rect(colour = "transparent"),
    axis.text=element_text(size=15))
p2swFungi

# Panel D: fungal cooling.
p2scFungi <- autoplot.prcWithoutSPBact(fit_fungiCooling, xlab = "Growing season time point", ylab = "Cooling effect on the\nabundance of fungi") +
  scale_colour_manual(values = c("navyblue", "lightblue", "navy", "lightblue", "navy")) +
  scale_linetype_manual(values = c("dashed", "dashed", "dotted", "solid", "solid")) +
  # scale_y_continuous(limits = c(-0.5, 0.2), breaks = seq(-0.4,0, by=0.2)) +
  theme(
    legend.position = c(0.2, 0.12), 
    legend.title = element_blank(), 
    panel.background=element_blank(), 
    axis.line = element_line(colour = "black"), 
    legend.key.width=unit(4,"line"), 
    legend.key = element_rect(colour = "transparent"),
    axis.text=element_text(size=15))
p2scFungi

# Combine fungal panels C and D for interactive inspection.
RDA_speciesfungi <- plot_grid(p2swFungi, p2scFungi, align= c("hv"), labels=c("A", "B"), label_size = 20)
RDA_speciesfungi


# Combine and save final Figure 3 panels A–D.
RDA_speciesmicrobe <- plot_grid(p2swBact, p2scBact, p2swFungi, p2scFungi, align= c("hv"), labels=c("A", "B", "C", "D"), label_size = 20)
RDA_speciesmicrobe
ggsave("output/figures/figure_3_microbial_prc_warming_cooling.png", width = 17, height= 17)
