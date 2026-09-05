# Manuscript table preparation
#
# Converts complete analysis outputs into presentation-ready tables used in the
# manuscript. Analysis tables remain unchanged in output/tables; formatted
# manuscript tables are written to output/tables/manuscript.

manuscript_table_dir <- "output/tables/manuscript"
dir.create(manuscript_table_dir, recursive = TRUE, showWarnings = FALSE)


# Table 2: Mean PRC treatment weights ------------------------------------

treatment_scores <- read.csv(
  "output/tables/cross_community_prc_treatment_score_summary.csv",
  check.names = FALSE
)

required_score_columns <- c(
  "community", "Treatment", "mean", "se", "significant"
)
if (!all(required_score_columns %in% names(treatment_scores))) {
  stop("The treatment-score input is missing columns required for Table 2.")
}

community_names <- c(
  cooled_bacteria_gradient = "Bacteria/Archaea",
  warmed_bacteria_gradient = "Bacteria/Archaea",
  cooled_fungi_gradient = "Fungi",
  warmed_fungi_gradient = "Fungi",
  cooled_plants_gradient = "Plants",
  warmed_plants_gradient = "Plants"
)

treatment_scores$Community <- unname(
  community_names[treatment_scores$community]
)
if (anyNA(treatment_scores$Community)) {
  stop("Table 2 contains an unrecognized community label.")
}

# An asterisk indicates that the absolute mean exceeds two standard errors,
# following the rule stated in the manuscript.
treatment_scores$formatted_mean <- sprintf("%.2f", treatment_scores$mean)
treatment_scores$formatted_mean[treatment_scores$significant] <- paste0(
  treatment_scores$formatted_mean[treatment_scores$significant],
  "*"
)
treatment_scores$formatted_se <- sprintf("%.2f", treatment_scores$se)

treatment_order <- c(
  "Low C2", "Low C1", "Mid C1", "High W2", "High W1", "Mid W1"
)
community_order <- c("Bacteria/Archaea", "Fungi", "Plants")

if (nrow(treatment_scores) != length(treatment_order) * length(community_order)) {
  stop("Table 2 should contain one row for each treatment-community combination.")
}

table_2 <- data.frame(Treatment = treatment_order, check.names = FALSE)
for (community_name in community_order) {
  community_scores <- treatment_scores[
    treatment_scores$Community == community_name,
    c("Treatment", "formatted_mean", "formatted_se")
  ]
  community_scores <- community_scores[
    match(treatment_order, community_scores$Treatment),
  ]

  if (anyNA(community_scores$Treatment)) {
    stop(paste("Table 2 is missing one or more treatments for", community_name))
  }

  table_2[[paste(community_name, "treatment weight")]] <-
    community_scores$formatted_mean
  table_2[[paste(community_name, "standard error")]] <-
    community_scores$formatted_se
}

write.csv(
  table_2,
  file.path(manuscript_table_dir, "table_2_prc_treatment_weights.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Table 3A: PRC taxon-weight correlations -------------------------------

weight_correlations <- read.csv(
  "output/tables/cross_community_prc_taxon_weight_correlations.csv",
  check.names = FALSE
)

required_correlation_columns <- c(
  "community", "direction", "correlation"
)
if (!all(required_correlation_columns %in% names(weight_correlations))) {
  stop("The correlation input is missing columns required for Table 3A.")
}

table_3_communities <- c("Plants", "Fungi", "Bacteria/Archaea")
if (nrow(weight_correlations) != 2 * length(table_3_communities)) {
  stop("Table 3A should contain warming and cooling results for each community.")
}

table_3a <- data.frame(Community = table_3_communities, check.names = FALSE)
for (direction_name in c("Warming", "Cooling")) {
  direction_results <- weight_correlations[
    weight_correlations$direction == direction_name,
    c("community", "correlation")
  ]
  direction_results <- direction_results[
    match(table_3_communities, direction_results$community),
  ]

  if (anyNA(direction_results$community)) {
    stop(paste("Table 3A is missing one or more", direction_name, "results."))
  }

  table_3a[[paste(direction_name, "vs. gradient, r")]] <- sprintf(
    "%.2f", direction_results$correlation
  )
}

write.csv(
  table_3a,
  file.path(manuscript_table_dir, "table_3a_prc_taxon_weight_correlations.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Table 3B: Microbial ASV overlap ---------------------------------------

asv_overlap <- read.csv(
  "output/tables/microbial_differential_abundance_gradient_overlap.csv",
  check.names = FALSE
)

required_overlap_columns <- c(
  "community", "direction", "experimental_asvs", "overlapping_asvs",
  "overlap_percent"
)
if (!all(required_overlap_columns %in% names(asv_overlap))) {
  stop("The ASV-overlap input is missing columns required for Table 3B.")
}

microbial_communities <- c("Bacteria/Archaea", "Fungi")
if (nrow(asv_overlap) != 2 * length(microbial_communities)) {
  stop("Table 3B should contain warming and cooling results for each microbial community.")
}

table_3b <- data.frame(Community = microbial_communities, check.names = FALSE)
for (direction_name in c("Warming", "Cooling")) {
  direction_results <- asv_overlap[
    asv_overlap$direction == direction_name,
    c("community", "experimental_asvs", "overlapping_asvs", "overlap_percent")
  ]
  direction_results <- direction_results[
    match(microbial_communities, direction_results$community),
  ]

  if (anyNA(direction_results$community)) {
    stop(paste("Table 3B is missing one or more", direction_name, "results."))
  }

  table_3b[[direction_name]] <- sprintf(
    "%d/%d (%d%%)",
    direction_results$overlapping_asvs,
    direction_results$experimental_asvs,
    round(direction_results$overlap_percent)
  )
}

write.csv(
  table_3b,
  file.path(manuscript_table_dir, "table_3b_microbial_asv_overlap.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Table S4: Plant cover and PRC taxon weights ---------------------------

plant_cover <- read.csv(
  "output/tables/plant_mean_scaled_cover_by_control.csv",
  check.names = FALSE
)
plant_weights <- read.csv(
  "output/tables/plant_prc_species_weight_correlations.csv",
  check.names = FALSE
)

required_cover_columns <- c(
  "species", "family", "High Control", "Mid Control", "Low Control"
)
required_weight_columns <- c(
  "species", "family", "warmed_plants_gradient", "warmed_plants",
  "cooled_plants_gradient", "cooled_plants"
)
if (!all(required_cover_columns %in% names(plant_cover))) {
  stop("The plant-cover input is missing columns required for Table S4.")
}
if (!all(required_weight_columns %in% names(plant_weights))) {
  stop("The plant-weight input is missing columns required for Table S4.")
}

plant_table_data <- merge(
  plant_weights,
  plant_cover[, setdiff(required_cover_columns, "family")],
  by = "species",
  all.x = TRUE,
  sort = FALSE
)

missing_families <- plant_table_data$species[is.na(plant_table_data$family)]
if (!identical(missing_families, "Unknown round leaves")) {
  stop("Unexpected missing or unmatched plant-family assignments in Table S4.")
}

# Table S4 includes identified taxa only, matching the manuscript.
plant_table_data <- plant_table_data[!is.na(plant_table_data$family), ]
plant_table_data <- plant_table_data[
  order(plant_table_data$warmed_plants_gradient, decreasing = TRUE),
]

format_cover <- function(x) {
  x[is.na(x)] <- 0
  sprintf("%.1f", 100 * x)
}

table_s4 <- data.frame(
  `Plant species` = plant_table_data$species,
  `Plant family` = plant_table_data$family,
  `High control mean proportional cover (%)` = format_cover(
    plant_table_data[["High Control"]]
  ),
  `Mid control mean proportional cover (%)` = format_cover(
    plant_table_data[["Mid Control"]]
  ),
  `Low control mean proportional cover (%)` = format_cover(
    plant_table_data[["Low Control"]]
  ),
  `Warmed gradient PRC weight` = sprintf(
    "%.2f", plant_table_data$warmed_plants_gradient
  ),
  `Warmed experiment PRC weight` = sprintf(
    "%.2f", plant_table_data$warmed_plants
  ),
  `Cooled gradient PRC weight` = sprintf(
    "%.2f", plant_table_data$cooled_plants_gradient
  ),
  `Cooled experiment PRC weight` = sprintf(
    "%.2f", plant_table_data$cooled_plants
  ),
  check.names = FALSE
)

if (nrow(table_s4) != 69) {
  stop("Table S4 should contain 69 identified plant taxa.")
}

write.csv(
  table_s4,
  file.path(manuscript_table_dir, "table_s4_plant_cover_and_prc_weights.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Table S5: Microbial PRC weights and control abundances ----------------

bacteria_weights <- read.csv(
  "output/tables/bacteria_archaea_prc_species_scores_with_taxonomy.csv",
  check.names = FALSE
)
fungi_weights <- read.csv(
  "output/tables/fungi_prc_species_weight_correlations.csv",
  check.names = FALSE
)

format_microbial_weights <- function(data, community, prefix) {
  required_columns <- c(
    "ASV_id", paste0("warmed_", prefix, "_gradient"),
    paste0("warmed_", prefix), paste0("cooled_", prefix, "_gradient"),
    paste0("cooled_", prefix), "High Control", "Mid Control", "Low Control",
    "Genus", "Family", "Phylum", "Class", "Order"
  )
  if (!all(required_columns %in% names(data))) {
    stop(paste(community, "data are missing columns required for Table S5."))
  }
  if (anyDuplicated(data$ASV_id)) {
    stop(paste(community, "data contain duplicate ASV identifiers."))
  }

  data.frame(
    community = community,
    ASV_id = data$ASV_id,
    LogAbundance_HighControl = round(data[["High Control"]], 2),
    LogAbundance_MidControl = round(data[["Mid Control"]], 2),
    LogAbundance_LowControl = round(data[["Low Control"]], 2),
    weight_warmed_gradient = round(
      data[[paste0("warmed_", prefix, "_gradient")]], 2
    ),
    weight_warmed_experiment = round(data[[paste0("warmed_", prefix)]], 2),
    weight_cooled_gradient = round(
      data[[paste0("cooled_", prefix, "_gradient")]], 2
    ),
    weight_cooled_experiment = round(data[[paste0("cooled_", prefix)]], 2),
    Genus = data$Genus,
    Family = data$Family,
    Phylum = data$Phylum,
    Class = data$Class,
    Order = data$Order,
    check.names = FALSE
  )
}

table_s5 <- rbind(
  format_microbial_weights(
    bacteria_weights, "bacteria and archaea", "bacteria"
  ),
  format_microbial_weights(fungi_weights, "fungi", "fungi")
)

write.csv(
  table_s5,
  file.path(manuscript_table_dir, "table_s5_microbial_prc_weights.csv"),
  row.names = FALSE,
  na = "NA",
  quote = FALSE
)


# Table S6: Highest absolute PRC responders -----------------------------
# Within each community, select the 15 largest absolute weights from each
# PRC category and retain the union. Ties are resolved deterministically by
# taxon identifier after ordering by absolute weight.

select_top_prc_responders <- function(data, id_column, weight_columns, n = 15) {
  selection_categories <- setNames(
    vector("list", nrow(data)),
    as.character(data[[id_column]])
  )

  for (category_name in names(weight_columns)) {
    weight_column <- weight_columns[[category_name]]
    eligible <- which(!is.na(data[[weight_column]]))
    ordered <- eligible[order(
      -abs(data[[weight_column]][eligible]),
      as.character(data[[id_column]][eligible])
    )]
    selected <- head(ordered, n)

    for (row_number in selected) {
      taxon_id <- as.character(data[[id_column]][row_number])
      selection_categories[[taxon_id]] <- c(
        selection_categories[[taxon_id]],
        category_name
      )
    }
  }

  selected_ids <- names(selection_categories)[lengths(selection_categories) > 0]
  selected_data <- data[match(selected_ids, as.character(data[[id_column]])), ]
  selected_data$selection_categories <- vapply(
    selection_categories[selected_ids],
    paste,
    collapse = "; ",
    FUN.VALUE = character(1)
  )
  selected_data
}

prc_categories <- c(
  "Warmed gradient" = "warmed_gradient",
  "Warmed experiment" = "warmed_experiment",
  "Cooled gradient" = "cooled_gradient",
  "Cooled experiment" = "cooled_experiment"
)

plant_s6_source <- data.frame(
  taxon_id = plant_table_data$species,
  taxon_name = plant_table_data$species,
  phylum = NA_character_,
  family = plant_table_data$family,
  high_control = 100 * ifelse(
    is.na(plant_table_data[["High Control"]]), 0,
    plant_table_data[["High Control"]]
  ),
  mid_control = 100 * ifelse(
    is.na(plant_table_data[["Mid Control"]]), 0,
    plant_table_data[["Mid Control"]]
  ),
  low_control = 100 * ifelse(
    is.na(plant_table_data[["Low Control"]]), 0,
    plant_table_data[["Low Control"]]
  ),
  warmed_gradient = plant_table_data$warmed_plants_gradient,
  warmed_experiment = plant_table_data$warmed_plants,
  cooled_gradient = plant_table_data$cooled_plants_gradient,
  cooled_experiment = plant_table_data$cooled_plants,
  stringsAsFactors = FALSE
)

make_microbial_s6_source <- function(data, prefix) {
  data.frame(
    taxon_id = as.character(data$ASV_id),
    taxon_name = data$Genus,
    phylum = data$Phylum,
    family = data$Family,
    high_control = data[["High Control"]],
    mid_control = data[["Mid Control"]],
    low_control = data[["Low Control"]],
    warmed_gradient = data[[paste0("warmed_", prefix, "_gradient")]],
    warmed_experiment = data[[paste0("warmed_", prefix)]],
    cooled_gradient = data[[paste0("cooled_", prefix, "_gradient")]],
    cooled_experiment = data[[paste0("cooled_", prefix)]],
    stringsAsFactors = FALSE
  )
}

bacteria_s6_source <- make_microbial_s6_source(
  bacteria_weights, "bacteria"
)
fungi_s6_source <- make_microbial_s6_source(fungi_weights, "fungi")

plant_s6 <- select_top_prc_responders(
  plant_s6_source, "taxon_id", prc_categories, n = 15
)
bacteria_s6 <- select_top_prc_responders(
  bacteria_s6_source, "taxon_id", prc_categories, n = 15
)
fungi_s6 <- select_top_prc_responders(
  fungi_s6_source, "taxon_id", prc_categories, n = 15
)

format_s6_section <- function(data, community, control_measure) {
  maximum_weight <- apply(
    abs(data[unname(prc_categories)]),
    1,
    max,
    na.rm = TRUE
  )
  data <- data[order(-maximum_weight, data$taxon_id), ]

  data.frame(
    Community = community,
    `Taxon ID` = data$taxon_id,
    `Taxon name` = data$taxon_name,
    Phylum = data$phylum,
    Family = data$family,
    `Control measure` = control_measure,
    `High control` = round(data$high_control, 2),
    `Mid control` = round(data$mid_control, 2),
    `Low control` = round(data$low_control, 2),
    `Warmed gradient PRC weight` = round(data$warmed_gradient, 2),
    `Warmed experiment PRC weight` = round(data$warmed_experiment, 2),
    `Cooled gradient PRC weight` = round(data$cooled_gradient, 2),
    `Cooled experiment PRC weight` = round(data$cooled_experiment, 2),
    `Selected among top 15 for` = data$selection_categories,
    check.names = FALSE
  )
}

table_s6 <- rbind(
  format_s6_section(
    plant_s6, "Plants", "Mean scaled proportional dominance (%)"
  ),
  format_s6_section(
    bacteria_s6, "Bacteria/Archaea", "Mean log-transformed abundance"
  ),
  format_s6_section(
    fungi_s6, "Fungi", "Mean log-transformed abundance"
  )
)

write.csv(
  table_s6,
  file.path(manuscript_table_dir, "table_s6_highest_prc_responders.csv"),
  row.names = FALSE,
  na = "NA",
  quote = FALSE
)


# Tables 1 and S3: Pairwise PRC tests ----------------------------------

pairwise_results <- read.csv(
  "output/tables/cross_community_pairwise_prc_results.csv",
  check.names = FALSE
)

required_pairwise_columns <- c(
  "community", "treatmentOriginGroupCompare", "eig RDA1", "Df", "F",
  "p value"
)
if (!all(required_pairwise_columns %in% names(pairwise_results))) {
  stop("The pairwise PRC input is missing columns required for Table S3.")
}
if (nrow(pairwise_results) != 45) {
  stop("Table S3 should contain 45 pairwise PRC results.")
}


# Table 1: Summary of community responses ------------------------------
# P values are retrieved directly from the complete pairwise PRC results.
# Response interpretations are retained explicitly because they also depend
# on the PRC trajectories shown in Figures 2 and 3, not only on significance.

table_1_community_names <- c(
  plants = "Plants",
  `bacteria and archaea` = "Bacteria/Archaea",
  fungi = "Fungi"
)
table_1_treatments <- c(
  "High W2", "High W1", "Mid W1", "Low C2", "Low C1", "Mid C1"
)

table_1_responses <- c(
  "Plants|High W2" =
    "Shifted toward destination; remained distinct from destination",
  "Plants|High W1" =
    "Shifted toward destination; remained distinct from destination",
  "Plants|Mid W1" =
    "Shifted toward destination; remained distinct from destination",
  "Plants|Low C2" =
    "Weak directional movement; remained statistically similar to origin",
  "Plants|Low C1" =
    "Weak directional movement; remained statistically similar to origin",
  "Plants|Mid C1" =
    "Persistent offset from both origin and destination",
  "Bacteria/Archaea|High W2" = "Converged toward destination",
  "Bacteria/Archaea|High W1" =
    "Unresolved; origin and destination controls did not differ",
  "Bacteria/Archaea|Mid W1" =
    "Directional movement toward destination; endpoint comparisons unresolved",
  "Bacteria/Archaea|Low C2" = "Remained close to origin",
  "Bacteria/Archaea|Low C1" =
    "Directional movement toward destination; endpoint comparisons unresolved",
  "Bacteria/Archaea|Mid C1" =
    "Unresolved; origin and destination controls did not differ",
  "Fungi|High W2" = "Weak, seasonally variable movement toward destination",
  "Fungi|High W1" = "Weak, seasonally variable movement toward destination",
  "Fungi|Mid W1" = "Weak, seasonally variable movement toward destination",
  "Fungi|Low C2" = "Remained similar to origin",
  "Fungi|Low C1" = "Remained similar to origin",
  "Fungi|Mid C1" =
    "Directional movement toward destination; endpoint comparisons unresolved"
)

get_table_1_p_value <- function(data, treatment, comparison_type) {
  pattern <- paste0(
    "^", treatment, " and .* \\(", comparison_type, "\\)$"
  )
  matching_rows <- grepl(pattern, data$treatmentOriginGroupCompare)

  if (sum(matching_rows) != 1) {
    stop(paste(
      "Table 1 expected one", comparison_type, "comparison for", treatment
    ))
  }

  data[["p value"]][matching_rows]
}

format_table_1_p_value <- function(value) {
  ifelse(value < 0.01, sprintf("%.3f", value), sprintf("%.2f", value))
}

table_1_rows <- list()
row_number <- 1

for (analysis_community in names(table_1_community_names)) {
  manuscript_community <- table_1_community_names[[analysis_community]]
  community_results <- pairwise_results[
    pairwise_results$community == analysis_community,
  ]

  for (treatment in table_1_treatments) {
    response_key <- paste(manuscript_community, treatment, sep = "|")
    table_1_rows[[row_number]] <- data.frame(
      Community = manuscript_community,
      Treatment = treatment,
      `Origin control P` = format_table_1_p_value(
        get_table_1_p_value(community_results, treatment, "origin")
      ),
      `Destination control P` = format_table_1_p_value(
        get_table_1_p_value(community_results, treatment, "destination")
      ),
      `PRC response` = unname(table_1_responses[response_key]),
      check.names = FALSE
    )
    row_number <- row_number + 1
  }
}

table_1 <- do.call(rbind, table_1_rows)
if (nrow(table_1) != 18 || anyNA(table_1)) {
  stop("Table 1 should contain 18 complete treatment-community rows.")
}

write.csv(
  table_1,
  file.path(manuscript_table_dir, "table_1_community_response_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Table S3: Pairwise PRC tests ------------------------------------------

table_s3 <- data.frame(
  Community = pairwise_results$community,
  `Treatment Comparison` = pairwise_results$treatmentOriginGroupCompare,
  `eig RDA1` = sprintf("%.2f", pairwise_results[["eig RDA1"]]),
  df = pairwise_results$Df,
  F = sprintf("%.2f", pairwise_results$F),
  `p value` = ifelse(
    pairwise_results[["p value"]] < 0.01,
    sprintf("%.3f", pairwise_results[["p value"]]),
    sprintf("%.2f", pairwise_results[["p value"]])
  ),
  check.names = FALSE
)

write.csv(
  table_s3,
  file.path(manuscript_table_dir, "table_s3_pairwise_prc_tests.csv"),
  row.names = FALSE,
  quote = FALSE
)
