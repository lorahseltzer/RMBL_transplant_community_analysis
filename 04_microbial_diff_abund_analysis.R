# Microbial differential-abundance analysis: RMBL transplant experiment
#
# 1. Bacteria/archaea (16S) differential abundance
# 2. Fungal (ITS) differential abundance
# 3. Figure 4 heatmap

library(DESeq2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(grid)

# Fix the random-number state for reproducible analysis and plotting behavior.
set.seed(20250904)

# Send any incidental graphics produced during a batch run to a null device.
# Named figure files are still written by their explicit ggsave() calls.
batch_graphics_device <- !interactive()
if (batch_graphics_device) {
  grDevices::pdf(file = NULL)
}

# Bacteria/archaea differential abundance -------------------------------

# Load filtered 16S data
topbact <- read.csv("data/bacteria_archaea_asv_table.csv")

# Create an unnormalized ASV-by-sample count matrix; DESeq2 estimates library
# size factors internally.
bactcount <- select(topbact, Sample, ASV_id, Abundance) %>% pivot_wider(names_from = Sample, values_from = Abundance,values_fill=0 ) %>% column_to_rownames("ASV_id")

# Create sample metadata and use treatment as a factor.
bactcolumn <- select(topbact, Sample, treatmentOrigin, season) %>% distinct(Sample, .keep_all = TRUE) %>% column_to_rownames("Sample")
bactcolumn$treatmentOrigin <- as.factor(sub(" ", "_", bactcolumn$treatmentOrigin)) # remove spaces
bactcolumn$season <- as.factor(bactcolumn$season)

# Align count-matrix columns with metadata rows.
bactcount <- bactcount[, rownames(bactcolumn)]
all(rownames(bactcolumn) == colnames(bactcount))

# Fit a DESeq2 model that accounts for season before testing treatment.
deseqbact <- DESeqDataSetFromMatrix(bactcount, bactcolumn, design=~season + treatmentOrigin)

# Use the default Wald significance test.
deseqbact <- DESeq(deseqbact)#, test = LRT, reduced=~treatmentOrigin )
res <- results(deseqbact)
res
mcols(res)$description # what the columns are.You can see that they're not reporting the results from the comparison I'm interested in.
resdf <- as.data.frame(res)

# Bacteria/archaea treatment contrasts -----------------------------------

### WARMED
### Which taxa have changed when warmed in comparison to the origin site? (High W2 vs High Control) ###
Res_deseqbact <- results(deseqbact, contrast=c("treatmentOrigin", "High_W2", "High_Control"), alpha=0.05)
# Log2 fold change = effect size, the ratio between the abundance in the contrast condition and its abundance in the reference condition. This shows the log2 fold change in High W2 relative to High Control. So if it is negative, it means the taxa decreased in High W2 relative to the High Control.
summary(Res_deseqbact)
Res_deseqbact
sum(Res_deseqbact$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed (i.e. High W2 vs High Control)
res_HW2 <- Res_deseqbact %>% data.frame() %>% filter(padj < 0.05) 
res_HW2 <- res_HW2 %>% mutate(ASV_id = row.names(res_HW2)) %>% mutate(condition = "High W2 vs High Control")

### Which taxa have changed when warmed in comparison to the origin site? (High W1 vs High Control) ###
Res_deseqbact2 <- results(deseqbact, contrast=c("treatmentOrigin", "High_W1", "High_Control"), alpha=0.05)
sum(Res_deseqbact2$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed (i.e. High W1 vs High Control)
res_HW1 <- Res_deseqbact2 %>% data.frame() %>% filter(padj < 0.05) 
res_HW1 <- res_HW1 %>% mutate(ASV_id = row.names(res_HW1)) %>% mutate(condition = "High W1 vs High Control")

### Which taxa have changed when warmed in comparison to the origin site? (Mid W1 vs Mid Control) ###
Res_deseqbact3 <- results(deseqbact, contrast=c("treatmentOrigin", "Mid_W1", "Mid_Control"), alpha=0.05)
sum(Res_deseqbact3$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed (i.e. Mid W1 vs Mid Control)
res_MW1 <- Res_deseqbact3 %>% data.frame() %>% filter(padj < 0.05) 
res_MW1 <- res_MW1 %>% mutate(ASV_id = row.names(res_MW1)) %>% mutate(condition = "Mid W1 vs Mid Control")


# Join these 3 dataframes to see how many unique ASV were enriched or depleted across all warming treatments. 
warm_unique <- bind_rows(list(res_HW2, res_HW1, res_MW1))
warm_unique$ASV_id <- as.integer(warm_unique$ASV_id) # change ASV_id to class = integer
# How many unique ASVs are there in the data frame?
n_distinct(warm_unique$ASV_id) # 78 unique ASVs changing with all warmed treatments
warm_unique2 <- pivot_wider(warm_unique, id_cols="ASV_id", names_from="condition", values_from = "log2FoldChange") %>% rename(HW2_HC = "High W2 vs High Control", HW1_HC = "High W1 vs High Control", MW1_MC = "Mid W1 vs Mid Control")
warm_unique3 <- warm_unique2 %>% 
  mutate(across(.cols = c(2:4), .fns = function(x) ifelse(x > 0, 1, -1))) %>%
  group_by(ASV_id) %>% 
  mutate(mean= mean(c(HW2_HC, HW1_HC, MW1_MC), na.rm=TRUE)) 
sum(warm_unique3$mean == 1, na.rm = TRUE)
sum(warm_unique3$mean == -1, na.rm = TRUE)
# any other number means that ASV is enriched or depleted depending on treatment
# Count the number of enriched (1's), depleted (-1's), or both (anything other number) in this df


### CONTROLS: decreasing elevation (naturally warmed)
### Which taxa change naturally as elevation decreases along gradient? (Low Control vs High Control) ###
Res_deseqbact4 <- results(deseqbact, contrast=c("treatmentOrigin", "Low_Control", "High_Control"), alpha=0.05)
sum(Res_deseqbact4$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LCHC <- Res_deseqbact4 %>% data.frame() %>% filter(padj < 0.05) 
res_LCHC <- res_LCHC %>% mutate(ASV_id = row.names(res_LCHC)) %>% mutate(condition = "Low Control vs High Control")


### Which taxa change naturally as elevation decreases along gradient? (Low Control vs Mid Control) ###
Res_deseqbact5 <- results(deseqbact, contrast=c("treatmentOrigin", "Low_Control", "Mid_Control"), alpha=0.05)
sum(Res_deseqbact5$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LCMC <- Res_deseqbact5 %>% data.frame() %>% filter(padj < 0.05) 
res_LCMC <- res_LCMC %>% mutate(ASV_id = row.names(res_LCMC)) %>% mutate(condition = "Low Control vs Mid Control")

### Which taxa change naturally as elevation decreases along gradient? (Mid Control vs High Control) ###
Res_deseqbact6 <- results(deseqbact, contrast=c("treatmentOrigin", "Mid_Control", "High_Control"), alpha=0.05)
sum(Res_deseqbact6$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_MCHC <- Res_deseqbact6 %>% data.frame() %>% filter(padj < 0.05) 
res_MCHC <- res_MCHC %>% mutate(ASV_id = row.names(res_MCHC)) %>% mutate(condition = "Mid Control vs High Control")

### COOLED
### Which taxa have changed when cooled in comparison to the origin site? (Low C2 vs Low Control) ###
Res_deseqbact7 <- results(deseqbact, contrast=c("treatmentOrigin", "Low_C2", "Low_Control"), alpha=0.05)
sum(Res_deseqbact7$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LC2 <- Res_deseqbact7 %>% data.frame() %>% filter(padj < 0.05) 
res_LC2 <- res_LC2 %>% mutate(ASV_id = row.names(res_LC2)) %>% mutate(condition = "Low C2 vs Low Control")

### Which taxa have changed when cooled in comparison to the origin site? (Low C1 vs Low Control) ###
Res_deseqbact8 <- results(deseqbact, contrast=c("treatmentOrigin", "Low_C1", "Low_Control"), alpha=0.05)
sum(Res_deseqbact8$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LC1 <- Res_deseqbact8 %>% data.frame() %>% filter(padj < 0.05) 
res_LC1 <- res_LC1 %>% mutate(ASV_id = row.names(res_LC1)) %>% mutate(condition = "Low C1 vs Low Control")

### Which taxa have changed when cooled in comparison to the origin site? (Mid C1 vs Mid Control) ###
Res_deseqbact9 <- results(deseqbact, contrast=c("treatmentOrigin", "Mid_C1", "Mid_Control"), alpha=0.05)
sum(Res_deseqbact9$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_MC1 <- Res_deseqbact9 %>% data.frame() %>% filter(padj < 0.05) 
res_MC1 <- res_MC1 %>% mutate(ASV_id = row.names(res_MC1)) %>% mutate(condition = "Mid C1 vs Mid Control")

# Join these 3 dataframes to see how many unique ASV were enriched or depleted across all cooling treatments. 
cool_unique <- bind_rows(list(res_LC2, res_LC1, res_MC1))
cool_unique$ASV_id <- as.integer(cool_unique$ASV_id) # change ASV_id to class = integer
# How many unique ASVs are there in the data frame?
n_distinct(cool_unique$ASV_id) # 67 unique ASVs changing with all cooled treatments
cool_unique2 <- pivot_wider(cool_unique, id_cols="ASV_id", names_from="condition", values_from = "log2FoldChange") %>% rename(LC2_LC = "Low C2 vs Low Control", LC1_LC = "Low C1 vs Low Control", MC1_MC = "Mid C1 vs Mid Control")
cool_unique3 <- cool_unique2 %>% 
  mutate(across(.cols = c(2:4), .fns = function(x) ifelse(x > 0, 1, -1))) %>%
  group_by(ASV_id) %>% 
  mutate(mean= mean(c(LC2_LC, LC1_LC, MC1_MC), na.rm=TRUE)) 
sum(cool_unique3$mean == 1, na.rm = TRUE)
sum(cool_unique3$mean == -1, na.rm = TRUE)
# Count the number of enriched (1's), depleted (-1's), or both (anything other number) in this df


### CONTROLS: increasing elevation (naturally cooled)
### Which taxa change naturally as elevation increases along gradient? (High Control vs Low Control) ###
Res_deseqbact10 <- results(deseqbact, contrast=c("treatmentOrigin", "High_Control", "Low_Control"), alpha=0.05)
sum(Res_deseqbact10$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LCHCC <- Res_deseqbact10 %>% data.frame() %>% filter(padj < 0.05) 
res_LCHCC <- res_LCHCC %>% mutate(ASV_id = row.names(res_LCHCC)) %>% mutate(condition = "High Control vs Low Control")

### Which taxa change naturally as elevation increases along gradient? (Mid Control vs Low Control) ###
Res_deseqbact11 <- results(deseqbact, contrast=c("treatmentOrigin", "Mid_Control", "Low_Control"), alpha=0.05)
sum(Res_deseqbact11$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LCMCC <- Res_deseqbact11 %>% data.frame() %>% filter(padj < 0.05) 
res_LCMCC <- res_LCMCC %>% mutate(ASV_id = row.names(res_LCMCC)) %>% mutate(condition = "Mid Control vs Low Control")

### Which taxa change naturally as elevation increases along gradient? (High Control vs Mid Control) ###
Res_deseqbact12 <- results(deseqbact, contrast=c("treatmentOrigin", "High_Control", "Mid_Control"), alpha=0.05)
sum(Res_deseqbact12$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_MCHCC <- Res_deseqbact12 %>% data.frame() %>% filter(padj < 0.05) 
res_MCHCC <- res_MCHCC %>% mutate(ASV_id = row.names(res_MCHCC)) %>% mutate(condition = "High Control vs Mid Control")

# another way to look at results is with the name command. #One exception to the equivalence of these two commands, is that, using contrast will additionally set to 0 the estimated LFC in a comparison of two groups, where all of the counts in the two groups are equal to 0 (while other groups have positive counts). As this may be a desired feature to have the LFC in these cases set to 0, one can use contrast to build these results tables
resultsNames(deseqbact) # lists the coefficients
#res <- results(deseqbact, name="treatmentOrigin_High_W2_vs_High_Control")

# Combine all of the data frames containing differentially abundant ASVs that you created above
res_diffabund <- bind_rows(list(res_HW2, res_HW1, res_MW1, res_LCHC, res_LCMC, res_MCHC, res_LC2, res_LC1, res_MC1, res_LCHCC, res_LCMCC, res_MCHCC))
res_diffabund$ASV_id <- as.integer(res_diffabund$ASV_id) # change ASV_id to class = integer
# How many unique ASVs are there in the data frame?
n_distinct(res_diffabund$ASV_id)

# Make a dataframe of ASV_id's for each OTU sequence from topbact (also taxonomic info for future use)
topbactASV <- topbact %>% distinct(topbact$ASV_id, .keep_all=TRUE) %>% select("ASV_id", "OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
# Join the OTU sequence to the differential abundance dataframe based on ASV_id, remove duplicates
res_diffabund2 <- inner_join(res_diffabund, topbactASV, by="ASV_id") %>% distinct(ASV_id, condition, .keep_all=TRUE)
#write.csv(res_diffabund2, "output/tables/res_diffabund2.csv", row.names = FALSE)

# Filter by only the significantly abundant ASVs
abundASV <- filter(topbact, ASV_id %in% res_diffabund2$ASV_id)
#write.csv(abundASV, "output/tables/abundASV.csv", row.names = FALSE)
#abundASV <- read.csv("output/tables/abundASV.csv")

# Show exploratory abundance histograms only during an interactive R session.
if (interactive()) {
  hist(
    topbact$relAbundTreat,
    main = "Average relative abundance of all ASVs by treatment",
    xlab = "Relative abundance"
  )
  hist(topbact$absAbundTreat, xlab = "Absolute Abundance")
  hist(
    abundASV$relAbundTreat,
    main = "Relative abundance of 128 ASVs changing significantly with treatments",
    xlab = "Relative abundance"
  )
}

### Add abundance by treatment columns to differential abundance results (res_diffabund2): the representation of each ASV in each of the pairwise compared treatments

# Make columns for condition 1, condition 2 of the pairwise comparison by parsing column "condition". Add keys for join by ASV_id and condition.
res_diffabund2 <- res_diffabund2 %>% 
  separate(condition, c("condition1", "condition2"), sep=" vs ", remove=FALSE) %>% 
  unite(ASV_condition1, c("ASV_id", "condition1"), remove = FALSE) %>% 
  unite(ASV_condition2, c("ASV_id", "condition2"), remove = FALSE) 

# Remove all but relevant columns in abundASV. Make keys for join. Remove duplicates.
abundASV_reduce <- abundASV %>% 
  select("ASV_id", "treatmentOrigin", "relAbundTreat", "absAbundTreat") %>% 
  unite(ASV_condition1, c("ASV_id", "treatmentOrigin")) %>% 
  mutate(ASV_condition2 = ASV_condition1) %>% 
  distinct(ASV_condition1, .keep_all = TRUE)

# Join abundASV$relAbundTreatment by ASV and condition 1 to make relAbundCondition1
res_diffabund3 <- left_join(res_diffabund2, abundASV_reduce, by="ASV_condition1") %>% select(-"ASV_condition2.y")
res_diffabund3 <- rename(res_diffabund3, "ASV_condition2" = "ASV_condition2.x")
res_diffabund3 <- rename(res_diffabund3, "relAbundCondition1" = "relAbundTreat")
res_diffabund3 <- rename(res_diffabund3, "absAbundCondition1" = "absAbundTreat")

# Join abundASV$relAbundTreatment by ASV and condition 2 to make relAbundCondition2
res_diffabund4 <- left_join(res_diffabund3, abundASV_reduce, by = "ASV_condition2")%>% select(-"ASV_condition1.y") 
res_diffabund4 <- rename(res_diffabund4, "ASV_condition1" = "ASV_condition1.x")
res_diffabund4 <- rename(res_diffabund4, "relAbundCondition2" = "relAbundTreat")
res_diffabund4 <- rename(res_diffabund4, "absAbundCondition2" = "absAbundTreat")

# Change NA to zero in relAbundCondition1 and relAbundCondition2 to show that the ASV wasn't present.
res_diffabund4$relAbundCondition1[is.na(res_diffabund4$relAbundCondition1)] <- 0
res_diffabund4$relAbundCondition2[is.na(res_diffabund4$relAbundCondition2)] <- 0
res_diffabund4$absAbundCondition1[is.na(res_diffabund4$absAbundCondition1)] <- 0
res_diffabund4$absAbundCondition2[is.na(res_diffabund4$absAbundCondition2)] <- 0

# Add a "key" column that combines ASV_id, Phylum, and Family for use with heatmap below.
res_diffabund4 <- res_diffabund4 %>% 
  unite(key, ASV_id, Phylum, sep = "   ", remove=FALSE) %>% 
  unite(key1, key, Family, sep = ", ", remove=FALSE) %>% 
  unite(key2, ASV_id, Family, sep = "  ")
write.csv(
  res_diffabund4,
  "output/tables/bacteria_archaea_differential_abundance_results.csv",
  row.names = FALSE
)



# Fungal differential abundance -----------------------------------------

# Load filtered ITS data created by "LS_1_FilterPhyloseqData_ITS.R"
topfungi <- read.csv("data/fungal_asv_table.csv")


# Create an unnormalized ASV-by-sample count matrix; DESeq2 estimates library
# size factors internally.
fungicount <- select(topfungi, Sample, ASV_id, Abundance) %>% pivot_wider(names_from = Sample, values_from = Abundance,values_fill=0 ) %>% column_to_rownames("ASV_id")

# Create sample metadata and use treatment and season as factors.
fungicolumn <- select(topfungi, Sample, treatmentOrigin, season) %>% distinct(Sample, .keep_all = TRUE) %>% column_to_rownames("Sample")
fungicolumn$treatmentOrigin <- as.factor(sub(" ", "_", fungicolumn$treatmentOrigin)) # remove spaces
fungicolumn$season <- as.factor(fungicolumn$season)

# Align count-matrix columns with metadata rows.
fungicount <- fungicount[, rownames(fungicolumn)]
all(rownames(fungicolumn) == colnames(fungicount))

# Fit a DESeq2 model that accounts for season before testing treatment.
deseqfungi <- DESeqDataSetFromMatrix(fungicount, fungicolumn, design=~season + treatmentOrigin)

# Use the default Wald significance test.
deseqfungi <- DESeq(deseqfungi)#, test = LRT, reduced=~treatmentOrigin )
resf <- results(deseqfungi)
resf
mcols(resf)$description # what the columns are.You can see that they're not reporting the results from the comparison I'm interested in
resdff <- as.data.frame(resf)

# The deseq2 function uses the Benjamini-Hochberg method to calculate the adjusted p-value by default. There are options to use other methods in the results function. But we want to use the q-value/Storey method which is not an option. You need to use the separate qvalue function. *** I did this with the below contrasts and got the same number of signficant ASVs as with using the adjusted p-value <0.05, so decided to skip this. ***
#browseVignettes("qvalue")
#pvalues <- res$pvalue #create a list of pvalues from your results
#qres <- qvalue(p=pvalues)
#summary(qres) # There are 122(bact) taxa with p-value <0.025. This means that there is a 2.5% chance of false positives, so with 1663 spots * 0.025, we expect 41.6 false positives. There are 50 taxa with a q-value <0.025. This means that out of those 50, we should expect 50*0.025 = 1.25 false positives. There are 55 taxa with a q-value <0.05, 55*0.05=2.75 false positives. Better to use the cut off at 0.025.
#hist(qres)
#plot(qres)
#qvalues <- as.data.frame(qres$qvalues) # make a dataframe of just the q-values if you want to look closer




# Fungal treatment contrasts --------------------------------------------

# With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor, the comparison will be the last level of this variable over the reference level. However, the order of the variables of the design do not matter so long as the user specifies the comparison to build a results table for, using the name or contrast arguments of results.

### WARMED
### Which taxa have changed when warmed in comparison to the origin site? (High W2 vs High Control) ###
Res_deseqfungi <- results(deseqfungi, contrast=c("treatmentOrigin", "High_W2", "High_Control"), alpha=0.05)
# Log2 fold change = effect size, the ratio between the abundance in the contrast condition and its abundance in the reference condition. This shows the log2 fold change in High W2 relative to High Control. So if it is negative, it means the taxa decreased in High W2 relative to the High Control.
summary(Res_deseqfungi)
Res_deseqfungi
sum(Res_deseqfungi$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed (i.e. High W2 vs High Control)
res_HW2f <- Res_deseqfungi %>% data.frame() %>% filter(padj < 0.05) 
res_HW2f <- res_HW2f %>% mutate(ASV_id = row.names(res_HW2f)) %>% mutate(condition = "High W2 vs High Control")

### Which taxa have changed when warmed in comparison to the origin site? (High W1 vs High Control) ###
Res_deseqfungi2 <- results(deseqfungi, contrast=c("treatmentOrigin", "High_W1", "High_Control"), alpha=0.05)
sum(Res_deseqfungi2$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed (i.e. High W1 vs High Control)
res_HW1f <- Res_deseqfungi2 %>% data.frame() %>% filter(padj < 0.05) 
res_HW1f <- res_HW1f %>% mutate(ASV_id = row.names(res_HW1f)) %>% mutate(condition = "High W1 vs High Control")

### Which taxa have changed when warmed in comparison to the origin site? (Mid W1 vs Mid Control) ###
Res_deseqfungi3 <- results(deseqfungi, contrast=c("treatmentOrigin", "Mid_W1", "Mid_Control"), alpha=0.05)
sum(Res_deseqfungi3$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed (i.e. Mid W1 vs Mid Control)
res_MW1f <- Res_deseqfungi3 %>% data.frame() %>% filter(padj < 0.05) 
res_MW1f <- res_MW1f %>% mutate(ASV_id = row.names(res_MW1f)) %>% mutate(condition = "Mid W1 vs Mid Control")


# Join these 3 dataframes to see how many unique ASV were enriched or depleted across all warming treatments. 
warm_uniquef <- bind_rows(list(res_HW2f, res_HW1f, res_MW1f))
warm_uniquef$ASV_id <- as.integer(warm_uniquef$ASV_id) # change ASV_id to class = integer
# How many unique ASVs are there in the data frame?
n_distinct(warm_uniquef$ASV_id) # 67 unique ASVs changing with all cooled treatments
warm_uniquef2 <- pivot_wider(warm_uniquef, id_cols="ASV_id", names_from="condition", values_from = "log2FoldChange") %>% rename(HW2_HC = "High W2 vs High Control", HW1_HC = "High W1 vs High Control", MW1_MC = "Mid W1 vs Mid Control")
warm_uniquef3 <- warm_uniquef2 %>% 
  mutate(across(.cols = c(2:4), .fns = function(x) ifelse(x > 0, 1, -1))) %>%
  group_by(ASV_id) %>% 
  mutate(mean= mean(c(HW2_HC, HW1_HC, MW1_MC), na.rm=TRUE)) 
sum(warm_uniquef3$mean == 1, na.rm = TRUE)
sum(warm_uniquef3$mean == -1, na.rm = TRUE)
# Count the number of enriched (1's), depleted (-1's), or both (anything else) in this df



### CONTROLS: decreasing elevation (naturally warmed)
### Which taxa change naturally as elevation decreases along gradient? (Low Control vs High Control) ###
Res_deseqfungi4 <- results(deseqfungi, contrast=c("treatmentOrigin", "Low_Control", "High_Control"), alpha=0.05)
sum(Res_deseqfungi4$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LCHCf <- Res_deseqfungi4 %>% data.frame() %>% filter(padj < 0.05) 
res_LCHCf <- res_LCHCf %>% mutate(ASV_id = row.names(res_LCHCf)) %>% mutate(condition = "Low Control vs High Control")

### Which taxa change naturally as elevation decreases along gradient? (Low Control vs Mid Control) ###
Res_deseqfungi5 <- results(deseqfungi, contrast=c("treatmentOrigin", "Low_Control", "Mid_Control"), alpha=0.05)
sum(Res_deseqfungi5$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LCMCf <- Res_deseqfungi5 %>% data.frame() %>% filter(padj < 0.05) 
res_LCMCf <- res_LCMCf %>% mutate(ASV_id = row.names(res_LCMCf)) %>% mutate(condition = "Low Control vs Mid Control")

### Which taxa change naturally as elevation decreases along gradient? (Mid Control vs High Control) ###
Res_deseqfungi6 <- results(deseqfungi, contrast=c("treatmentOrigin", "Mid_Control", "High_Control"), alpha=0.05)
sum(Res_deseqfungi6$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_MCHCf <- Res_deseqfungi6 %>% data.frame() %>% filter(padj < 0.05) 
res_MCHCf <- res_MCHCf %>% mutate(ASV_id = row.names(res_MCHCf)) %>% mutate(condition = "Mid Control vs High Control")

### COOLED
### Which taxa have changed when cooled in comparison to the origin site? (Low C2 vs Low Control) ###
Res_deseqfungi7 <- results(deseqfungi, contrast=c("treatmentOrigin", "Low_C2", "Low_Control"), alpha=0.05)
sum(Res_deseqfungi7$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LC2f <- Res_deseqfungi7 %>% data.frame() %>% filter(padj < 0.05) 
res_LC2f <- res_LC2f %>% mutate(ASV_id = row.names(res_LC2f)) %>% mutate(condition = "Low C2 vs Low Control")

### Which taxa have changed when cooled in comparison to the origin site? (Low C1 vs Low Control) ###
Res_deseqfungi8 <- results(deseqfungi, contrast=c("treatmentOrigin", "Low_C1", "Low_Control"), alpha=0.05)
sum(Res_deseqfungi8$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LC1f <- Res_deseqfungi8 %>% data.frame() %>% filter(padj < 0.05) 
res_LC1f <- res_LC1f %>% mutate(ASV_id = row.names(res_LC1f)) %>% mutate(condition = "Low C1 vs Low Control")

### Which taxa have changed when cooled in comparison to the origin site? (Mid C1 vs Mid Control) ###
Res_deseqfungi9 <- results(deseqfungi, contrast=c("treatmentOrigin", "Mid_C1", "Mid_Control"), alpha=0.05)
sum(Res_deseqfungi9$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_MC1f <- Res_deseqfungi9 %>% data.frame() %>% filter(padj < 0.05) 
res_MC1f <- res_MC1f %>% mutate(ASV_id = row.names(res_MC1f)) %>% mutate(condition = "Mid C1 vs Mid Control")


# Join these 3 dataframes to see how many unique ASV were enriched or depleted across all cooling treatments. 
cool_uniquef <- bind_rows(list(res_LC2f, res_LC1f, res_MC1f))
cool_uniquef$ASV_id <- as.integer(cool_uniquef$ASV_id) # change ASV_id to class = integer
# How many unique ASVs are there in the data frame?
n_distinct(cool_uniquef$ASV_id) # 146 unique ASVs changing with all cooled treatments
cool_unique2f <- pivot_wider(cool_uniquef, id_cols="ASV_id", names_from="condition", values_from = "log2FoldChange") %>% rename(LC2_LC = "Low C2 vs Low Control", LC1_LC = "Low C1 vs Low Control", MC1_MC = "Mid C1 vs Mid Control")
cool_unique3f <- cool_unique2f %>% 
  mutate(across(.cols = c(2:4), .fns = function(x) ifelse(x > 0, 1, -1))) %>%
  group_by(ASV_id) %>% 
  mutate(mean= mean(c(LC2_LC, LC1_LC, MC1_MC), na.rm=TRUE)) 
sum(cool_unique3f$mean == 1, na.rm = TRUE)
sum(cool_unique3f$mean == -1, na.rm = TRUE)
# Count the number of enriched (1's), depleted (-1's), or both (anything other number) in this df

### CONTROLS: increasing elevation (naturally cooled)
### Which taxa change naturally as elevation increases along gradient? (High Control vs Low Control) ###
Res_deseqfungi10 <- results(deseqfungi, contrast=c("treatmentOrigin", "High_Control", "Low_Control"), alpha=0.05)
sum(Res_deseqfungi10$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LCHCCf <- Res_deseqfungi10 %>% data.frame() %>% filter(padj < 0.05) 
res_LCHCCf <- res_LCHCCf %>% mutate(ASV_id = row.names(res_LCHCCf)) %>% mutate(condition = "High Control vs Low Control")

### Which taxa change naturally as elevation increases along gradient? (Mid Control vs Low Control) ###
Res_deseqfungi11 <- results(deseqfungi, contrast=c("treatmentOrigin", "Mid_Control", "Low_Control"), alpha=0.05)
sum(Res_deseqfungi11$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_LCMCCf <- Res_deseqfungi11 %>% data.frame() %>% filter(padj < 0.05) 
res_LCMCCf <- res_LCMCCf %>% mutate(ASV_id = row.names(res_LCMCCf)) %>% mutate(condition = "Mid Control vs Low Control")

### Which taxa change naturally as elevation increases along gradient? (High Control vs Mid Control) ###
Res_deseqfungi12 <- results(deseqfungi, contrast=c("treatmentOrigin", "High_Control", "Mid_Control"), alpha=0.05)
sum(Res_deseqfungi12$padj < 0.05, na.rm=TRUE) # how many ASVs have adjusted p < 0.05 (Benjamini-Hochberg)
# Make data frame of differentially abundant ASV_IDs and a column labeling when they changed
res_MCHCCf <- Res_deseqfungi12 %>% data.frame() %>% filter(padj < 0.05) 
res_MCHCCf <- res_MCHCCf %>% mutate(ASV_id = row.names(res_MCHCCf)) %>% mutate(condition = "High Control vs Mid Control")

# another way to look at results is with the name command. #One exception to the equivalence of these two commands, is that, using contrast will additionally set to 0 the estimated LFC in a comparison of two groups, where all of the counts in the two groups are equal to 0 (while other groups have positive counts). As this may be a desired feature to have the LFC in these cases set to 0, one can use contrast to build these results tables
resultsNames(deseqfungi) # lists the coefficients
#res <- results(deseqbact, name="treatmentOrigin_High_W2_vs_High_Control")

# Combine all of the data frames containing differentially abundant ASVs that you created above
res_diffabundf <- bind_rows(list(res_HW2f, res_HW1f, res_MW1f, res_LCHCf, res_LCMCf, res_MCHCf, res_LC2f, res_LC1f, res_MC1f, res_LCHCCf, res_LCMCCf, res_MCHCCf))
res_diffabundf$ASV_id <- as.integer(res_diffabundf$ASV_id) # change ASV_id to class = integer
# How many unique ASVs are there in the data frame?
n_distinct(res_diffabundf$ASV_id)


# Experimental and elevation-gradient ASV overlap -----------------------
# Count significant ASVs from the experimental warming or cooling
# contrasts that also changed significantly across the corresponding
# natural elevation contrasts. These counts are summarized in Table 3B.

summarize_asv_overlap <- function(
    community, direction, experimental_results, gradient_results) {
  experimental_asvs <- unique(experimental_results$ASV_id)
  gradient_asvs <- unique(gradient_results$ASV_id)
  overlapping_asvs <- intersect(experimental_asvs, gradient_asvs)

  data.frame(
    community = community,
    direction = direction,
    experimental_asvs = length(experimental_asvs),
    overlapping_asvs = length(overlapping_asvs),
    overlap_percent = 100 * length(overlapping_asvs) /
      length(experimental_asvs)
  )
}

microbial_asv_overlap <- bind_rows(
  summarize_asv_overlap(
    "Bacteria/Archaea", "Warming",
    bind_rows(res_HW2, res_HW1, res_MW1),
    bind_rows(res_LCHC, res_LCMC, res_MCHC)
  ),
  summarize_asv_overlap(
    "Bacteria/Archaea", "Cooling",
    bind_rows(res_LC2, res_LC1, res_MC1),
    bind_rows(res_LCHCC, res_LCMCC, res_MCHCC)
  ),
  summarize_asv_overlap(
    "Fungi", "Warming",
    bind_rows(res_HW2f, res_HW1f, res_MW1f),
    bind_rows(res_LCHCf, res_LCMCf, res_MCHCf)
  ),
  summarize_asv_overlap(
    "Fungi", "Cooling",
    bind_rows(res_LC2f, res_LC1f, res_MC1f),
    bind_rows(res_LCHCCf, res_LCMCCf, res_MCHCCf)
  )
)

write.csv(
  microbial_asv_overlap,
  "output/tables/microbial_differential_abundance_gradient_overlap.csv",
  row.names = FALSE
)

# Make a dataframe of ASV_id's for each OTU sequence from topbact (also taxonomic info for future use)
topfungiASV <- topfungi %>% distinct(topfungi$ASV_id, .keep_all=TRUE) %>% select("ASV_id", "OTU", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Join the OTU sequence to the differential abundance dataframe based on ASV_id, remove duplicates
res_diffabund2f <- inner_join(res_diffabundf, topfungiASV, by="ASV_id") %>% distinct(ASV_id, condition, .keep_all=TRUE)
#write.csv(res_diffabund2f, "output/tables/res_diffabund2f_ITS.csv", row.names = FALSE)

#######################
# CALCULATE ABUNDANCES AT TREATMENT LEVEL #

# Goal: How representative of the microbial community of each treatment are the differentially abundant ASVs? 
# I have relativeabundances by sample,  but need relative abundance by treatment because I compared differentially abundant taxa by treatment. I calculated this in "LS_1_FilterPhyloseqData_ITS.R" in the column "relAbundTreat" aka topfungi$relAbundTreat

# Filter by only the significantly abundant ASVs
abundASVf <- filter(topfungi, ASV_id %in% res_diffabund2f$ASV_id)
#write.csv(abundASV, "output/tables/abundASV.csv", row.names = FALSE)

# Show exploratory abundance histograms only during an interactive R session.
if (interactive()) {
  hist(
    topfungi$relAbundTreat,
    main = "Average relative abundance of all ASVs by treatment",
    xlab = "Relative abundance"
  )
  hist(
    abundASVf$relAbundTreat,
    main = "Relative abundance of 285 ASVs changing significantly with treatments",
    xlab = "Relative abundance"
  )
}

### Add abundance by treatment columns to differential abundance results (res_diffabund2): the representation of each ASV in each of the pairwise compared treatments

# Make columns for condition 1, condition 2 of the pairwise comparison by parsing column "condition". Add keys for join = ASV_id and condition.
res_diffabund2f <- res_diffabund2f %>% 
  separate(condition, c("condition1", "condition2"), sep=" vs ", remove=FALSE) %>% 
  unite(ASV_condition1, c("ASV_id", "condition1"), remove = FALSE) %>% 
  unite(ASV_condition2, c("ASV_id", "condition2"), remove = FALSE) 

# Remove all but relevant columns in abundASV. Make keys for join. Remove duplicates.
abundASV_reducef <- abundASVf %>% select("ASV_id", "treatmentOrigin", "relAbundTreat") %>% 
  unite(ASV_condition1, c("ASV_id", "treatmentOrigin")) %>% 
  mutate(ASV_condition2 = ASV_condition1) %>% # a place holder to make column?
  distinct(ASV_condition1, .keep_all = TRUE) #remove duplicates

# Join abundASV$relAbundTreatment by ASV and condition 1 to make relAbundCondition1
res_diffabund3f <- left_join(res_diffabund2f, abundASV_reducef, by="ASV_condition1") %>% 
  select(-"ASV_condition2.y") #remove this column since redundant
res_diffabund3f <- dplyr::rename(res_diffabund3f, "ASV_condition2"= "ASV_condition2.x")
res_diffabund3f <- dplyr::rename(res_diffabund3f, "relAbundCondition1" = "relAbundTreat")

# Join abundASV$relAbundTreatment by ASV and condition 2 to make relAbundCondition2
res_diffabund4f <- left_join(res_diffabund3f, abundASV_reducef, by = "ASV_condition2") %>% 
  select(-"ASV_condition1.y")
res_diffabund4f <- dplyr::rename(res_diffabund4f, "ASV_condition1" = "ASV_condition1.x")
res_diffabund4f <- dplyr::rename(res_diffabund4f, "relAbundCondition2" = "relAbundTreat")

# Change NA to zero in relAbundCondition1 and relAbundCondition2 to show that the ASV wasn't present.
res_diffabund4f$relAbundCondition1[is.na(res_diffabund4f$relAbundCondition1)] <- 0
res_diffabund4f$relAbundCondition2[is.na(res_diffabund4f$relAbundCondition2)] <- 0

# Add a "key" column that combines ASV_id, Phylum, and Family for use with heatmap below.
res_diffabund4f <- res_diffabund4f %>% 
  unite(key, ASV_id, Phylum, sep = "   ", remove=FALSE) %>% 
  unite(key1, key, Family, sep = ", ", remove=FALSE) %>% 
  unite(key2, ASV_id, Family, sep = "  ", remove=FALSE)

# Make a column "species" with the Genus and specific epithet (for use in merging by species when assigning traits)
res_diffabund4f <- res_diffabund4f %>%
  unite(species, Genus, Species, sep = "_", remove=FALSE) 
res_diffabund4f$species <- na_if(res_diffabund4f$species, "NA_NA")
write.csv(
  res_diffabund4f,
  "output/tables/fungi_differential_abundance_results.csv",
  row.names = FALSE
)



# Figure 4: differential-abundance heatmaps ------------------------------
# Bacteria/archaea filtering and heatmap.
### Filter the 128 ASVs by those that are most abundant (>0.1% of community in either condition). Then filter by log2fold change >10 or <-10.  

res_diffabund5 <- res_diffabund4 %>% filter(relAbundCondition1 > 0.001 | relAbundCondition2 > 0.001) # leaves 92 unique ASVs
res_diffabund6 <- res_diffabund5 %>% filter(log2FoldChange > 10 | log2FoldChange < -10) # leaves 63 unique ASVs


### Make a heatmap of log2 fold change by treatment ###

# Make a matrix of log2fold changes with row names as ASVs and treatment contrasts as columns.  !
# !!!! Can change which res_diffabund based on how you want the data filtered in the heatmap !!!!!
# This matrix reflects control contrasts that are naturally warmed down elevation.
log2df <- res_diffabund6 %>% select(key2, log2FoldChange, condition) %>% spread(condition, log2FoldChange) %>%column_to_rownames("key2") %>% replace(is.na(.), 0)
#log2df <- log2df[, c("Low Control vs High Control", "High Control vs Low Control", "High W2 vs High Control", "Low C2 vs Low Control", "Low Control vs Mid Control", "Mid Control vs Low Control", "Mid W1 vs Mid Control", "Low C1 vs Low Control", "Mid Control vs High Control", "High Control vs Mid Control", "High W1 vs High Control", "Mid C1 vs Mid Control")] # includes "naturally cooled" order of controls
log2dfnew <- log2df[, c("Low Control vs High Control", "High W2 vs High Control", "Low C2 vs Low Control", "Low Control vs Mid Control", "Mid W1 vs Mid Control", "Low C1 vs Low Control", "Mid Control vs High Control", "High W1 vs High Control", "Mid C1 vs Mid Control")]
log2dfnew2<- log2dfnew %>% dplyr::rename("Low vs High Control" = "Low Control vs High Control", "High W2 vs Origin" = "High W2 vs High Control", "Low C2 vs Origin" = "Low C2 vs Low Control", "Low vs Mid Control" = "Low Control vs Mid Control", "Mid W1 vs Origin" = "Mid W1 vs Mid Control", "Low C1 vs Origin" = "Low C1 vs Low Control", "Mid vs High Control" = "Mid Control vs High Control", "High W1 vs Origin" = "High W1 vs High Control", "Mid C1 vs Origin" = "Mid C1 vs Mid Control") %>% filter_all(any_vars(. != 0)) # Get names into shape for figures, remove rows where all are equal to zero
log2ma <- data.matrix(log2dfnew2)
class(log2ma)


contrastnew <- c("Low Control vs High Control", 
                 "High W2 vs High Control", 
                 "Low C2 vs Low Control", 
                 "Low Control vs Mid Control", 
                 "Mid W1 vs Mid Control", 
                 "Low C1 vs Low Control", 
                 "Mid Control vs High Control", 
                 "High W1 vs High Control", 
                 "Mid C1 vs Mid Control")
contrastnewsimple <- c("Low vs High Control", 
                       "High W2 vs Origin",
                       "Low C2 vs Origin", 
                       "Low vs Mid Control", 
                       "Mid W1 vs Origin",
                       "Low C1 vs Origin", 
                       "Mid vs High Control", 
                       "High W1 vs Origin", 
                       "Mid C1 vs Origin")
splitnew <- c("naturally warmed", 
              "transplant warmed",
              "transplant cooled", 
              "naturally warmed", 
              "transplant warmed",
              "transplant cooled",
              "naturally warmed",
              "transplant warmed",
              "transplant cooled")
comparison <- c(1,1,1,2,2,2,3,3,3)
splitdatanew <- data.frame(contrastnew,contrastnewsimple, splitnew, comparison)

# Third, use complexHeatmap to make the plot
# in order of treatments
htnew <- Heatmap(log2ma, name = "log2 fold change", column_title = "treatment contrast", column_title_side = "bottom", row_title="ASVs", cluster_columns=FALSE, column_split = splitdatanew$comparison, column_gap = grid::unit(5, "mm"))
if (interactive()) {
  draw(htnew)
}



### ITS Fungi

### Filter the 285 ASVs by those that are most abundant (>0.5% of community in either condition (for 16S we did >0.1% but this results in a large number with ITS)). Then filter by log2fold change >10 or <-10.  

res_diffabund5f <- res_diffabund4f %>% filter(relAbundCondition1 > 0.005 | relAbundCondition2 > 0.005) 
n_distinct(res_diffabund5f$OTU) # leaves 109 unique ASVs
res_diffabund6f <- res_diffabund5f %>% filter(log2FoldChange > 10 | log2FoldChange < -10) 
n_distinct(res_diffabund6f$OTU) # leaves 64 unique ASVs

# Show exploratory filtering plots only during an interactive R session.
if (interactive()) {
  hist(
    res_diffabund4f$log2FoldChange,
    xlab = "log2 fold change",
    main = "Log2 fold change - all 128 ASVs"
  )
  hist(
    res_diffabund5f$log2FoldChange,
    xlab = "log2 fold change",
    main = "Log2 fold change - 109 most abundant ASVs"
  )
  hist(
    res_diffabund6f$log2FoldChange,
    xlab = "log2 fold change",
    main = "Log2 fold change - 64 most abundant ASVs with largest log2"
  )
  plot(res_diffabund5f$log2FoldChange, res_diffabund5f$relAbundCondition1)
  plot(res_diffabund5f$relAbundCondition1, res_diffabund5f$relAbundCondition2)
  plot(res_diffabund6f$relAbundCondition1, res_diffabund6f$relAbundCondition2)
}

write.csv(
  res_diffabund6,
  "output/tables/bacteria_archaea_figure_4_filtered_asvs.csv",
  row.names = FALSE
)

# Fungal filtering and heatmap.

# Make a matrix of log2fold changes with row names as ASVs and treatment contrasts as columns.  !
# !!!! Can change which res_diffabund based on how you want the data filtered in the heatmap !!!!!
# This matrix reflects control contrasts that are naturally warmed down elevation.
log2dff <- res_diffabund6f %>% 
  select(key2, log2FoldChange, condition) %>% 
  spread(condition, log2FoldChange) %>%
  column_to_rownames("key2") %>% 
  replace(is.na(.), 0)
log2dffnew <- log2dff[, c("Low Control vs High Control", "High W2 vs High Control", "Low C2 vs Low Control", "Low Control vs Mid Control", "Mid W1 vs Mid Control", "Low C1 vs Low Control", "Mid Control vs High Control", "High W1 vs High Control", "Mid C1 vs Mid Control")]
log2dffnew2<- log2dffnew %>% dplyr::rename("Low vs High Control" = "Low Control vs High Control", "High W2 vs Origin" = "High W2 vs High Control", "Low C2 vs Origin" = "Low C2 vs Low Control", "Low vs Mid Control" = "Low Control vs Mid Control", "Mid W1 vs Origin" = "Mid W1 vs Mid Control", "Low C1 vs Origin" = "Low C1 vs Low Control", "Mid vs High Control" = "Mid Control vs High Control", "High W1 vs Origin" = "High W1 vs High Control", "Mid C1 vs Origin" = "Mid C1 vs Mid Control") %>% filter_all(any_vars(. != 0)) # Get names into shape for figures, remove rows where all are equal to zero
log2maf <- data.matrix(log2dffnew2)
class(log2maf)

contrastnew <- c("Low Control vs High Control", 
                 "High W2 vs High Control", 
                 "Low C2 vs Low Control", 
                 "Low Control vs Mid Control", 
                 "Mid W1 vs Mid Control", 
                 "Low C1 vs Low Control", 
                 "Mid Control vs High Control", 
                 "High W1 vs High Control", 
                 "Mid C1 vs Mid Control")
contrastnewsimple <- c("Low vs High Control", 
                       "High W2 vs Origin",
                       "Low C2 vs Origin", 
                       "Low vs Mid Control", 
                       "Mid W1 vs Origin",
                       "Low C1 vs Origin", 
                       "Mid vs High Control", 
                       "High W1 vs Origin", 
                       "Mid C1 vs Origin")
splitnew <- c("naturally warmed", 
              "transplant warmed",
              "transplant cooled", 
              "naturally warmed", 
              "transplant warmed",
              "transplant cooled",
              "naturally warmed",
              "transplant warmed",
              "transplant cooled")
comparison <- c(1,1,1,2,2,2,3,3,3)
splitdatanew <- data.frame(contrastnew,contrastnewsimple, splitnew, comparison)

# Third, use complexHeatmap to make the plot
# in order of treatments
htnewf <- Heatmap(log2maf, name = "log2 fold change", column_title = "treatment contrast", column_title_side = "bottom", row_title="ASVs", cluster_columns=FALSE, column_split = splitdatanew$comparison, column_gap = grid::unit(5, "mm"))
if (interactive()) {
  draw(htnewf)
}

write.csv(
  res_diffabund6f,
  "output/tables/fungi_figure_4_filtered_asvs.csv",
  row.names = FALSE
)


# Combine bacteria/archaea and fungal heatmaps into Figure 4.

# Convert the ComplexHeatmap objects to grobs for assembly with cowplot.
# Padding order is bottom, left, top, and right.
htnew2 <- grid.grabExpr(
  draw(htnew, padding = grid::unit(c(2, 2, 30, 2), "mm"))
)
htnewf2 <- grid.grabExpr(
  draw(htnewf, padding = grid::unit(c(2, 2, 30, 2), "mm"))
)

bacteria_heatmap_panel <- ggdraw(htnew2) +
  draw_label("Bacteria/Archaea", x = 0.5, y = 0.99, vjust = 1,
             fontface = "bold", size = 30)

fungi_heatmap_panel <- ggdraw(htnewf2) +
  draw_label("Fungi", x = 0.5, y = 0.99, vjust = 1,
             fontface = "bold", size = 30)

figure_4 <- plot_grid(
  bacteria_heatmap_panel,
  fungi_heatmap_panel,
  labels = c("A", "B"),
  label_size = 30,
  ncol = 1
)
if (interactive()) {
  print(figure_4)
}

ggsave(
  "output/figures/figure_4_microbial_differential_abundance_heatmap.png",
  plot = figure_4,
  width = 13,
  height = 24,
  units = "in",
  dpi = 300
)

if (batch_graphics_device) {
  grDevices::dev.off()
}
