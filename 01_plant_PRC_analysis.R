# Plant community analysis: RMBL transplant experiment
#
# 1. Plant Principal Response Curve (PRC) analysis
# 2. Correlations among plant PRC species weights
# 3. Figure 2 and associated species-weight linestacks

library(vegan)
library(devtools)
library(ggvegan)
library(ggplot2)
library(precrec)
library(tidyverse)
library(cowplot)
library(broom)

# Make permutation tests reproducible.
set.seed(20250901)

# Data preparation --------------------------------------------------------

# Get 2017 plant abundance data
abund2017 <- read.csv(
  "data/plant_abundance_2017.csv",
  header = TRUE,
  na.strings = "NA"
)
abund2017 <- abund2017 %>%
  filter(count != 0)

# Get plot metadata
plotInfo <- read.csv(
  "data/plant_plot_metadata.csv",
  header = TRUE,
  na.strings = "NA"
)

# Plant-family assignments transcribed from manuscript Table S4. The lookup
# preserves the family classification used in the manuscript; unknown plants
# are intentionally left without an inferred family.
plant_taxonomy <- read.csv(
  "data/plant_taxonomy.csv",
  header = TRUE,
  na.strings = "NA"
)
if (anyDuplicated(plant_taxonomy$species)) {
  stop("Plant taxonomy contains duplicate species names.")
}

# Abundance was measured in one turf per block. Expand those observations to
# all turfs in the corresponding origin-site block using the plot metadata.
abund2017expand <- abund2017 %>% 
  select(-turfID, -originPlotID) %>% # Remove turfID and originPlotID columns
  full_join(plotInfo, by = "originSiteBlock", relationship = "many-to-many") %>%
  filter(count != "NA") %>%
  select(-originSite.x, -treatment.x, -treatmentOrigin.x) %>%
  mutate(year = 2017)
abund2017expand <- dplyr::rename(
  abund2017expand,
  c(
    "originSite" = "originSite.y",
    "treatment" = "treatment.y",
    "treatmentOrigin" = "treatmentOrigin.y"
  )
)

# Get plant cover data
pctcover <- read.csv(
  "data/plant_community_cover_2018_2023.csv",
  header = TRUE,
  na.strings = NA
)
pctcover$percentCover[pctcover$percentCover == "F"] <- NA # Occurrence code, not percent cover
pctcover$percentCover <- as.integer(pctcover$percentCover)
coverall <- pctcover %>%
  filter(functionalGroup != "ground cover") %>%
  filter(
    turfID != "NA",
    !is.na(percentCover),
    !species %in% c("unknown seedling", "unknown forb", "unknown grass")
  ) %>%
  select(turfID, year, destinationSite, treatmentOriginGroup, functionalGroup, percentCover, species) %>%
  ungroup()
coverall$percentCover <- as.numeric(coverall$percentCover)

# Put 2017 abundance and 2018–2023 cover on a common proportional-dominance scale.

# Scale cover relative to the maximum cover within each turf and year.
coverscaled <- coverall %>%
  group_by(turfID, year) %>%
  mutate(scaled = percentCover / max(percentCover))

# Scale abundance relative to the maximum count within each turf.
abund2017scaled <- abund2017expand %>%
  group_by(turfID) %>%
  mutate(scaled = count / max(count)) %>%
  select(turfID, year, destinationSite, treatmentOriginGroup, scaled, count, species)

# Combine abundance and cover observations.
scaledcommunity <- bind_rows(coverscaled, abund2017scaled)

# Remove block 6 because it was transplanted one year later than the other blocks.
species_thin <- scaledcommunity %>% 
  filter(!turfID%in%c("mo6-1_c2_um6-3","mo6-2_c2_um6-5","mo6-3_c1_pf6-3","mo6-4_c1_pf6-4","mo6-5_c2_um6-6","pf6-1_c1_um6-2","pf6-2_c1_um6-1","pf6-4_w1_mo6-4","pf6-3_w1_mo6-3","um6-1_w1_pf6-2","um6-2_w1_pf6-1","um6-3_w2_mo6-1","um6-6_w2_mo6-5","um6-5_ws_um6-4","um6-4_w2_mo6-2")) %>%
  select(turfID, year, destinationSite, treatmentOriginGroup, species, scaled) %>%
  droplevels()

write.csv(species_thin, "output/tables/plant_community_scaled_long.csv", row.names = FALSE)

# Principal Response Curve analysis --------------------------------------


# Warming treatments. The first factor level is the control baseline.

speciesFat_Warming <- species_thin %>% 
  arrange(year) %>%
  filter(treatmentOriginGroup %in% c("High Control", "High W2", "High W1", "Mid W1", "Mid Control", "Low Control")) %>%
  spread(key = species, value = scaled, fill = 0) %>% 
  mutate(treatmentOriginGroup = as.character(treatmentOriginGroup)) %>% 
  mutate(year = factor(year), treatmentOriginGroup = factor(treatmentOriginGroup, levels = c("High Control", "High W1", "High W2", "Mid W1", "Mid Control", "Low Control"))) %>%
  ungroup()

# Retain only numeric community-response columns.
speciesData_Warming <- speciesFat_Warming %>%
  select(-turfID, -year, -destinationSite, -treatmentOriginGroup)

fit_speciesWarming <- prc(response = speciesData_Warming, treatment = speciesFat_Warming$treatmentOriginGroup, time = speciesFat_Warming$year, scale=TRUE)

summary(fit_speciesWarming)
plot(fit_speciesWarming, scaling="species", ylab="Warming effect on proportional \n dominance of plant species", xlab= "Year")



# Cooling treatments. The first factor level is the control baseline.

speciesFat_Cooling <- species_thin %>% 
  arrange(year) %>%
  filter(treatmentOriginGroup %in% c("Low Control", "Low C2", "Low C1", "Mid C1", "High Control", "Mid Control")) %>%
  spread(key = species, value = scaled, fill = 0) %>% 
  mutate(treatmentOriginGroup = as.character(treatmentOriginGroup)) %>% 
  mutate(year = factor(year), treatmentOriginGroup = factor(treatmentOriginGroup, levels = c("Low Control", "High Control", "Mid Control", "Mid C1","Low C1", "Low C2")))

# Retain only numeric community-response columns.
speciesData_Cooling <- speciesFat_Cooling %>%
  ungroup() %>%
  select(-turfID, -year, -destinationSite, -treatmentOriginGroup)

fit_speciesCooling <- prc(response = speciesData_Cooling, treatment = speciesFat_Cooling$treatmentOriginGroup, time = speciesFat_Cooling$year, scale = TRUE)

sumfitC <- summary(fit_speciesCooling)
plot(fit_speciesCooling, scaling="species", ylab="Cooling effect on proportional \n dominance of plant species", xlab= "Year")

# Extract species weights from the gradient PRCs before the correlation analyses
# that use them.
speciesscoresdf <- fortify(fit_speciesWarming) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community = "warmed_plants_gradient") %>%
  rownames_to_column(var = "species")

speciesscoresdfC <- fortify(fit_speciesCooling) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community = "cooled_plants_gradient") %>%
  rownames_to_column(var = "species")


# Correlations among plant PRC species weights ---------------------------

# Compare each species' experimental warming or cooling weight with its weight
# along the elevation gradient. Also compare warming and cooling weights to test
# for symmetry in species responses to temperature change.

# PRC species scores for the high-origin control and warming treatments.
HW2_PRC <- species_thin %>% 
  arrange(year) %>%
  filter(treatmentOriginGroup %in% c("High Control", "High W2", "High W1", "Mid W1")) %>%
  spread(key = species, value = scaled, fill = 0) %>% 
  mutate(treatmentOriginGroup = as.character(treatmentOriginGroup)) %>% 
  mutate(year = factor(year), treatmentOriginGroup = factor(treatmentOriginGroup, levels = c("High Control", "High W2", "High W1", "Mid W1"))) %>%
  ungroup()
# Make object containing only numerical cover data. 
Data_HW2_PRC<- HW2_PRC  %>%
  select(-turfID, -year, -destinationSite, -treatmentOriginGroup)
# Fit PRC.
fit_HW2_PRC  <- prc(response = Data_HW2_PRC, treatment = HW2_PRC$treatmentOriginGroup, time = HW2_PRC$year, scale=TRUE)
# Extract species scores.
scores_HW2_PRC<- fortify(fit_HW2_PRC) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="warmed_plants")%>%
  rownames_to_column(var="species")

# PRC species scores for the low-origin control and cooling treatments.
LC2_PRC <- species_thin %>% 
  arrange(year) %>%
  filter(treatmentOriginGroup %in% c("Low Control", "Low C2", "Low C1", "Mid C1")) %>%
  spread(key = species, value = scaled, fill = 0) %>% 
  mutate(treatmentOriginGroup = as.character(treatmentOriginGroup)) %>% 
  mutate(year = factor(year), treatmentOriginGroup = factor(treatmentOriginGroup, levels = c("Low Control", "Low C2", "Low C1", "Mid C1"))) %>%
  ungroup()
# Make object containing only numerical cover data. 
Data_LC2_PRC<- LC2_PRC  %>%
  select(-turfID, -year, -destinationSite, -treatmentOriginGroup)
# Fit PRC.
fit_LC2_PRC  <- prc(response = Data_LC2_PRC, treatment = LC2_PRC $treatmentOriginGroup, time = LC2_PRC$year, scale=TRUE)
# Extract species scores.
scores_LC2_PRC<- fortify(fit_LC2_PRC) %>%
  filter(score == "Species") %>%
  arrange(Response) %>%
  mutate(X = 1, community="cooled_plants")%>%
  rownames_to_column(var="species")

# Test for symmetry between warming and cooling responses.
weight_LC2_HW2 <- bind_rows(scores_LC2_PRC, scores_HW2_PRC)
weight_LC2_HW2 <- weight_LC2_HW2 %>% dplyr::rename(coverResponse = Response) %>% 
  pivot_wider(names_from = community, values_from = coverResponse)
cor.test(weight_LC2_HW2$cooled_plants, weight_LC2_HW2$warmed_plants)
plot(weight_LC2_HW2$cooled_plants, weight_LC2_HW2 $warmed_plants)

# Compare experimental warming with the natural warming gradient.
weight_HW_grad <- bind_rows(speciesscoresdf, scores_HW2_PRC)
weight_HW_grad <- weight_HW_grad %>% dplyr::rename(coverResponse = Response) %>%
  pivot_wider(names_from = community, values_from = coverResponse)
cor.test(weight_HW_grad$warmed_plants, weight_HW_grad$warmed_plants_gradient)
plot(weight_HW_grad$warmed_plants, weight_HW_grad$warmed_plants_gradient)

# Compare experimental cooling with the natural cooling gradient.
weight_LC_grad <- bind_rows(speciesscoresdfC, scores_LC2_PRC)
weight_LC_grad$Response <- weight_LC_grad$Response * (-1) # to reverse signs of weights to indicate flipped axis
weight_LC_grad <- weight_LC_grad %>% dplyr::rename(coverResponse = Response) %>%
  pivot_wider(names_from = community, values_from = coverResponse)
cor.test(weight_LC_grad$cooled_plants, weight_LC_grad$cooled_plants_gradient)
plot(weight_LC_grad$cooled_plants, weight_LC_grad$cooled_plants_gradient)

weight_warmcool_grad <- full_join(weight_HW_grad, weight_LC_grad, by="species") %>%
  select(species, warmed_plants_gradient, warmed_plants, cooled_plants_gradient, cooled_plants) %>%
  left_join(plant_taxonomy, by = "species") %>%
  relocate(family, .after = species)

write.csv(weight_warmcool_grad, "output/tables/plant_prc_species_weight_correlations.csv", row.names = FALSE)


# Descriptive summaries --------------------------------------------------
# Calculate mean scaled cover while accounting for turfs where a species was absent.
nturfs <- plotInfo %>%
  group_by(treatmentOriginGroup) %>%
  summarize(n_turfs = n_distinct(turfID)*6) # how many turfs per treatment origin group? multiplied by 6 years of data
scaled_nturfs <- left_join(scaledcommunity, nturfs, by="treatmentOriginGroup") # join n_turfs to scaled cover data for use when calculating means (to incorporate zeros)
meanscaledcover <- scaled_nturfs %>%
  filter(!is.na(scaled)) %>%
  group_by(year, treatmentOriginGroup, species) %>%
  reframe(mean = sum(scaled)/n_turfs) %>%
  distinct(year, treatmentOriginGroup, species, .keep_all=TRUE) %>%
  pivot_wider(names_from = "year", values_from="mean") # calculate mean scaled cover per year, treatment, species

meanscaledcover_treatment <- scaled_nturfs %>%
  filter(!is.na(scaled)) %>%
  group_by(treatmentOriginGroup, species) %>%
  reframe(mean = sum(scaled)/n_turfs) %>%
  distinct(treatmentOriginGroup, species, .keep_all=TRUE) %>%
  filter(treatmentOriginGroup %in% c("High Control", "Mid Control", "Low Control")) %>%
  mutate_if(is.numeric, round, 5) %>% 
  pivot_wider(names_from = "treatmentOriginGroup", values_from="mean") %>%
  left_join(plant_taxonomy, by = "species") %>%
  relocate(family, .after = species) # calculate mean scaled cover per treatment and species across all years
write.csv(meanscaledcover_treatment, file = "output/tables/plant_mean_scaled_cover_by_control.csv", row.names = FALSE)


# Figure 2 plotting functions --------------------------------------------

# PRC plot with a species-score rug.
autoplot.prc <- function(object, select, xlab, ylab,
                         title = NULL, subtitle = NULL, caption = NULL,
                         legend.position = "top", ...) {
  ## fortify the model object
  fobj <- fortify(object)
  
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
  plt <- plt + geom_hline(yintercept = 0, color = "grey")
  ## add species rug
  plt <- plt +
    geom_rug(data = spp,
             sides = "r",
             mapping = aes_string(group = NULL, x = NULL,
                                  colour = NULL, linetype = NULL))
  ## add the coefficients
  plt <- plt + geom_line() +
    theme(legend.position = legend.position) +
    scale_x_continuous(breaks = as.numeric(TimeLevs), minor_breaks = NULL)
  
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


# PRC plot without species scores.
autoplot.prcWithoutSP <- function(object, select, xlab, ylab,
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
    scale_x_continuous(breaks = as.numeric(TimeLevs), minor_breaks = NULL)
  
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

# Custom PRC plot with a species-score rug.
autoplot.prcCustom <- function(object, select, xlab, ylab,
                               title = NULL, subtitle = NULL, caption = NULL,
                               legend.position = "top", ...) {
  ## fortify the model object
  fobj <- fortify(object)
  
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
  ## add species rug
  plt <- plt +
    geom_rug(data = spp,
             sides = "r",
             mapping = aes_string(group = NULL, x = NULL,
                                  colour = NULL, linetype = NULL))
  ## add the coefficients
  plt <- plt + geom_line(size = 1.5) +
    theme(legend.position = legend.position,
          text = element_text(size=20),
          axis.text = element_text(size = 20)) +
    scale_x_continuous(breaks = as.numeric(TimeLevs), minor_breaks = NULL)
  
  ## add labels
  if(missing(xlab)) {
    xlab <- 'Time'
  }
  if(missing(ylab)) {
    ylab <- 'Treatment'
  }
  plt <- plt + labs(x = xlab, y = ylab, title = title, subtitle = subtitle,
                    caption = caption) #+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  ## return
  plt
}


# Figure 2 ---------------------------------------------------------------
# Warming treatments.
p2sw <- autoplot.prcWithoutSP(fit_speciesWarming, xlab="Year", ylab = "Warming effect on the \n proportional dominance of plant species") +
  scale_colour_manual(values = c("orange", "red", "red", "orange", "red")) +
  scale_linetype_manual(values = c("dashed", "dashed", "dotted", "solid", "solid")) +
  scale_y_continuous( breaks = pretty(fortify(fit_speciesWarming)$Response, n = 20)) +
  theme(legend.position = c(0.8, 0.7), legend.title = element_blank(), panel.background=element_blank(), axis.line = element_line(colour = "black"), legend.key.width=unit(4,"line"), legend.key = element_rect(colour = "transparent", fill = "white"))
p2sw

p3sw <- fortify(fit_speciesWarming) %>%
  filter(score == "Species") %>%
  mutate(X = 1) %>%
  ggplot(aes(x = label, y = Response, label = label)) +
  geom_text(aes(x = X), size=2) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks = pretty(fortify(fit_speciesWarming)$Response, n = 6), trans="reverse") +
  labs(x = "", y = "") +
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) 
p3sw

# Cooling treatments.
p2sc <- autoplot.prcWithoutSP(fit_speciesCooling, xlab = "Year", ylab = "Cooling effect on the \n proportional dominance of plant species") +
  scale_colour_manual(values = c("navyblue", "lightblue", "navyblue", "lightblue", "navyblue")) +
  scale_linetype_manual(values = c("solid", "solid", "dotted", "dashed", "dashed")) +
  scale_y_continuous(breaks = pretty(fortify(fit_speciesCooling)$Response, n = 20)) +
  theme(
    legend.position = c(0.8, 0.45), 
    legend.title = element_blank(), 
    panel.background=element_blank(), 
    axis.line = element_line(colour = "black"), 
    legend.key.width=unit(4,"line"), 
    legend.key = element_rect(colour = "transparent"))
p2sc

p3sc <- fortify(fit_speciesCooling) %>%
  filter(score == "species") %>%
  mutate(X = 1) %>%
  ggplot(aes(x = X, y = (Response), label = label)) +
  geom_text(aes(x = X), size=2) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks = pretty(fortify(fit_speciesCooling)$Response, n = 5), trans="reverse") +
  labs(x = "", y = "") +
  theme(panel.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
p3sc

TraitRDA_speciesCooling <- gridExtra::grid.arrange(p2sc,
                                                   layout_matrix = rbind(c(2,2,2,2,2,3,3)))


# Combine warming and cooling PRCs into Figure 2 panels A and B.
RDA_species <- plot_grid(p2sw, p2sc, align= c("hv"), labels=c("A", "B"), label_size = 20, nrow=2, ncol=1)
RDA_species
ggsave("output/figures/figure_2_plant_prc_warming_cooling.png", width = 8, height= 15)


# Figure 2 species weights -----------------------------------------------

# Save warmed and cooled species-weight linestacks for use with Figure 2
png("output/figures/figure_2_plant_prc_species_weights.png",
    width = 8, height = 15, units = "in", res = 300)
par(mfrow = c(2, 1))
linestack(speciesscoresdf$Response, labels = speciesscoresdf$label,
          cex = 0.4, axis = TRUE, at = 0, side = "right")
linestack(speciesscoresdfC$Response, labels = speciesscoresdfC$label,
          cex = 0.4, axis = TRUE, at = 0, side = "right")
dev.off()

# Combine warming and cooling species scores.
allspeciesscores <- full_join(speciesscoresdf, speciesscoresdfC, by="species") # combine species scores with warmed scores
write.csv(allspeciesscores, file="output/tables/plant_prc_warming_cooling_species_scores.csv", row.names = FALSE)

coverweights <- bind_rows(speciesscoresdf, speciesscoresdfC)
coverweights <- coverweights %>% dplyr::rename(coverResponse = Response)
