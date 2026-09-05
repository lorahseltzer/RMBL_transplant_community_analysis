# Cross-community Principal Response Curve comparison
#
# Compares PRC treatment scores and pairwise PRC tests among plants,
# bacteria/archaea, and fungi. The component analyses are run first so this
# script does not depend on objects left in an interactive R session.

source("01_plant_PRC_analysis.R")
source("02_microbial_PRC_analysis.R")

# Make cross-community permutation tests reproducible independently of the
# random-number state left by the component analyses.
set.seed(20250903)

# Taxon-weight correlations ---------------------------------------------
# Save the Pearson correlations between experimental PRC taxon weights and
# the corresponding natural elevation-gradient weights. These statistics are
# summarized for the manuscript in Table 3A.

summarize_weight_correlation <- function(
    community, direction, experimental_weights, gradient_weights) {
  complete_rows <- complete.cases(experimental_weights, gradient_weights)
  correlation_test <- cor.test(
    experimental_weights[complete_rows],
    gradient_weights[complete_rows],
    method = "pearson"
  )

  data.frame(
    community = community,
    direction = direction,
    n_taxa = sum(complete_rows),
    correlation = unname(correlation_test$estimate),
    confidence_interval_lower = correlation_test$conf.int[1],
    confidence_interval_upper = correlation_test$conf.int[2],
    statistic = unname(correlation_test$statistic),
    degrees_of_freedom = unname(correlation_test$parameter),
    p_value = correlation_test$p.value
  )
}

cross_community_weight_correlations <- bind_rows(
  summarize_weight_correlation(
    "Plants", "Warming", weight_HW_grad$warmed_plants,
    weight_HW_grad$warmed_plants_gradient
  ),
  summarize_weight_correlation(
    "Plants", "Cooling", weight_LC_grad$cooled_plants,
    weight_LC_grad$cooled_plants_gradient
  ),
  summarize_weight_correlation(
    "Fungi", "Warming", fweight_wgrad$warmed_fungi,
    fweight_wgrad$warmed_fungi_gradient
  ),
  summarize_weight_correlation(
    "Fungi", "Cooling", fweight_cgrad$cooled_fungi,
    fweight_cgrad$cooled_fungi_gradient
  ),
  summarize_weight_correlation(
    "Bacteria/Archaea", "Warming", weight_wgrad$warmed_bacteria,
    weight_wgrad$warmed_bacteria_gradient
  ),
  summarize_weight_correlation(
    "Bacteria/Archaea", "Cooling", weight_cgrad$cooled_bacteria,
    weight_cgrad$cooled_bacteria_gradient
  )
)

write.csv(
  cross_community_weight_correlations,
  "output/tables/cross_community_prc_taxon_weight_correlations.csv",
  row.names = FALSE
)


# Treatment-score comparison ---------------------------------------------
# Treatment scores describe the magnitude of community responses over time.
# PRCs that include destination controls are used because their trends match
# origin-only PRCs while producing more consistent score directions.

# Plant warming scores.
Tspeciesscoresdf <- fortify(fit_speciesWarming) %>%
  filter(score == "Sample") %>%
  arrange(Response) %>%
  mutate(X = 1, community = "warmed_plants_gradient")
# Plant cooling scores.
TspeciesscoresdfC <- fortify(fit_speciesCooling) %>%
  filter(score == "Sample") %>%
  arrange(Response) %>%
  mutate(X = 1, community = "cooled_plants_gradient")


# Bacteria/archaea warming scores.
Tbactscoresdf <- fortify(fit_bactWarming) %>%
  filter(score == "Sample") %>%
  arrange(Response) %>%
  mutate(X = 1, community = "warmed_bacteria_gradient")
# Bacteria/archaea cooling scores.
TbactscoresdfC <- fortify(fit_bactCooling) %>%
  filter(score == "Sample") %>%
  arrange(Response) %>%
  mutate(X = 1, community = "cooled_bacteria_gradient")

# Fungal warming scores.
Tfungiscoresdf <- fortify(fit_fungiWarming) %>%
  filter(score == "Sample") %>%
  arrange(Response) %>%
  mutate(X = 1, community = "warmed_fungi_gradient")
# Fungal cooling scores.
TfungiscoresdfC <- fortify(fit_fungiCooling) %>%
  filter(score == "Sample") %>%
  arrange(Response) %>%
  mutate(X = 1, community = "cooled_fungi_gradient")

# Combine treatment scores across communities and temperature directions.
alltreatmentscores <- bind_rows(
  TspeciesscoresdfC,
  Tspeciesscoresdf,
  Tfungiscoresdf,
  TfungiscoresdfC,
  Tbactscoresdf,
  TbactscoresdfC
)

# Summarize treatment scores across years or seasonal time points.
ats_stats <- alltreatmentscores %>%
  group_by(community, Treatment) %>%
  summarise(
    mean = mean(Response),
    sd = sd(Response),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(significant = abs(mean) > 2 * se) %>%
  filter(!Treatment %in% c("High Control", "Mid Control", "Low Control"))

write.csv(
  ats_stats,
  "output/tables/cross_community_prc_treatment_score_summary.csv",
  row.names = FALSE
)


# Pairwise PRC tests ------------------------------------------------------
# Define the treatment pairs used in the permutation tests.
treatmentPairs <- data.frame(pair = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15),
                             treatmentOriginGroup = c("High Control", "Low Control", 
                                                      "High Control", "Mid Control", 
                                                      "Mid Control", "Low Control",
                                                      "High W2", "High Control",
                                                      "High W2", "Low Control",
                                                      "High W1", "High Control",
                                                      "High W1", "Mid Control",
                                                      "Mid W1", "Mid Control",
                                                      "Mid W1", "Low Control",
                                                      "Low C2", "Low Control",
                                                      "Low C2", "High Control",
                                                      "Low C1", "Low Control",
                                                      "Low C1", "Mid Control",
                                                      "Mid C1", "Mid Control",
                                                      "Mid C1", "High Control"))

# Plant community.
results_prc <- data.frame()
results_prc_eigen <- data.frame()
for (x in unique(treatmentPairs$pair)) {
  treatmentPairs2 <- treatmentPairs %>% filter(pair == x)
  fat_data <- filter(species_thin, treatmentOriginGroup %in% treatmentPairs2$treatmentOriginGroup) %>% 
    spread(key = species, value = scaled, fill = 0) %>% 
    mutate(year = factor(year), treatmentOriginGroup = factor(treatmentOriginGroup, levels = c("High Control", "Mid Control", "Low Control", "High W2", "High W1", "Mid W1", "Low C2", "Low C1", "Mid C1"))) %>%
    ungroup()
  # Make object containing only numerical cover data. 
  speciesData <- fat_data %>%
    select(-turfID, -year, -destinationSite, -treatmentOriginGroup)
  set.seed(1575)
  prc_fitofone <- prc(response = speciesData, treatment = fat_data$treatmentOriginGroup, time = fat_data$year, scale = TRUE)
  testPRC <- anova.cca(prc_fitofone, step = 1000) %>% add_column(pair = x)
  results_prc <- rbind(results_prc, testPRC) %>% group_by(pair) %>% distinct()
  prc_eigen1 <- data.frame(prc_fitofone$CCA$eig) %>% add_column(pair = x)
  prc_eigen1$RDA <- rownames(prc_eigen1)
  results_prc_eigen <- rbind(results_prc_eigen, prc_eigen1) %>% distinct()
}

# Add descriptive labels to the pairwise test results.
treatmentPairsGrouped <- data.frame (pair = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), treatmentOriginGroupCompare = c("High Control and Low Control", "High Control and Mid Control", "Mid Control and Low Control", "High W2 and High Control (origin)", "High W2 and Low Control (destination)", "High W1 and High Control (origin)", "High W1 and Mid Control (destination)", "Mid W1 and Mid Control (origin)", "Mid W1 and Low Control (destination)", "Low C2 and Low Control (origin)", "Low C2 and High Control (destination)", "Low C1 and Low Control (origin)", "Low C1 and Mid Control (destination)", "Mid C1 and Mid Control (origin)", "Mid C1 and High Control (destination)"))

plant_prc_results <- results_prc_eigen %>% 
  filter(RDA == "RDA1") %>%
  left_join(results_prc, by = "pair") %>%
  left_join(treatmentPairsGrouped, by = "pair") %>%
  rename("eig RDA1" = "prc_fitofone.CCA.eig", "p value" = "Pr(>F)")  %>% 
  select("treatmentOriginGroupCompare", "eig RDA1", "Df", "F", "p value") %>% 
  mutate(community = "plants") %>% 
  distinct(treatmentOriginGroupCompare, .keep_all = TRUE)

write.csv(
  plant_prc_results,
  "output/tables/plant_pairwise_prc_results.csv",
  row.names = FALSE
)


# Bacteria/archaea community.
results_prcbact <- data.frame()
results_prc_eigenbact <- data.frame()
for (x in unique(treatmentPairs$pair)) {
  treatmentPairs2 <- treatmentPairs %>% filter(pair == x)
  fat_data <- filter(topbactASV, treatmentOrigin %in% treatmentPairs2$treatmentOriginGroup) %>% 
    spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
    mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels =  c("High Control", "Mid Control", "Low Control", "High W2", "High W1", "Mid W1", "Low C2", "Low C1", "Mid C1")))%>%
    ungroup()
  # Make object containing only numerical cover data. 
  speciesData <- fat_data %>%
    select(-turfID, -season, -treatmentOrigin, -seasonnum)
  set.seed(1575)
  prc_fitofone <- prc(response = speciesData, treatment = fat_data$treatmentOrigin, time = fat_data$seasonnum, scale = TRUE)
  testPRC <- anova.cca(prc_fitofone, step = 1000) %>% add_column(pair = x)
  results_prcbact <- rbind(results_prcbact, testPRC) %>% group_by(pair) %>% distinct()
  prc_eigen1 <- data.frame(prc_fitofone$CCA$eig) %>% add_column(pair = x)
  prc_eigen1$RDA <- rownames(prc_eigen1)
  results_prc_eigenbact <- rbind(results_prc_eigenbact, prc_eigen1) %>% distinct()
}

bact_prc_results <- results_prc_eigenbact %>% 
  filter(RDA == "RDA1") %>%
  left_join(results_prcbact, by = "pair") %>%
  left_join(treatmentPairsGrouped, by = "pair") %>%
  rename("eig RDA1" = "prc_fitofone.CCA.eig", "p value" = "Pr(>F)")  %>% 
  select("treatmentOriginGroupCompare", "eig RDA1", "Df", "F", "p value") %>% 
  mutate(community = "bacteria and archaea") %>% 
  distinct(treatmentOriginGroupCompare, .keep_all = TRUE)

write.csv(
  bact_prc_results,
  "output/tables/bacteria_archaea_pairwise_prc_results.csv",
  row.names = FALSE
)

# Fungal community.
results_prcfungi <- data.frame()
results_prc_eigenfungi <- data.frame()
for (x in unique(treatmentPairs$pair)) {
  treatmentPairs2 <- treatmentPairs %>% filter(pair == x)
  fat_data <- filter(topfungiASV, treatmentOrigin %in% treatmentPairs2$treatmentOriginGroup) %>% 
    spread(key = ASV_id, value = logAbundance, fill = 0) %>% 
    mutate(seasonnum = factor(seasonnum), treatmentOrigin = factor(treatmentOrigin, levels =  c("High Control", "Mid Control", "Low Control", "High W2", "High W1", "Mid W1", "Low C2", "Low C1", "Mid C1")))%>%
    ungroup()
  # Make object containing only numerical cover data. 
  speciesData <- fat_data %>%
    select(-turfID, -season, -treatmentOrigin, -seasonnum)
  set.seed(1575)
  prc_fitofone <- prc(response = speciesData, treatment = fat_data$treatmentOrigin, time = fat_data$seasonnum, scale = TRUE)
  testPRC <- anova.cca(prc_fitofone, step = 1000) %>% add_column(pair = x)
  results_prcfungi <- rbind(results_prcfungi, testPRC) %>% group_by(pair) %>% distinct()
  prc_eigen1 <- data.frame(prc_fitofone$CCA$eig) %>% add_column(pair = x)
  prc_eigen1$RDA <- rownames(prc_eigen1)
  results_prc_eigenfungi <- rbind(results_prc_eigenfungi, prc_eigen1) %>% distinct()
}

fungi_prc_results <- results_prc_eigenfungi %>% 
  filter(RDA == "RDA1") %>%
  left_join(results_prcfungi, by = "pair") %>%
  left_join(treatmentPairsGrouped, by = "pair") %>%
  rename("eig RDA1" = "prc_fitofone.CCA.eig", "p value" = "Pr(>F)")  %>% 
  select("treatmentOriginGroupCompare", "eig RDA1", "Df", "F", "p value") %>% 
  mutate(community = "fungi") %>% 
  distinct(treatmentOriginGroupCompare, .keep_all = TRUE)

write.csv(
  fungi_prc_results,
  "output/tables/fungi_pairwise_prc_results.csv",
  row.names = FALSE
)

# Combine pairwise PRC results across communities.
all_prc_results <- bind_rows(plant_prc_results, bact_prc_results, fungi_prc_results)
write.csv(
  all_prc_results,
  "output/tables/cross_community_pairwise_prc_results.csv",
  row.names = FALSE
)
