# Load necessary libraries for data manipulation, analysis, and visualization
library(mixOmics)
library(labdsv)
library(tidyverse)
library(Maaslin2)
library(BiocManager)
library(phyloseq)
library(ape)
library(forcats)
library(viridis)
library(knitr)
library(vegan)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ANCOMBC)
library(VennDiagram)

# Set working directory to the location of your data files
# Please update the path to your specific directory
setwd("/path/to/your/directory")

# Read in the data files: metabolomics data, taxonomy, and metadata
metabolomics <- as.data.frame(t(read.csv("metabolomics.csv", row.names = 1)))
taxonomy <- read.csv("taxon.csv", row.names = 1)
map <- read.csv("metadata.csv", row.names = 1)

# Quality control on metabolomics data
metabolomics.sums_fecal <- colSums(metabolomics)
# Filter metabolites with counts greater than 5
metabo_fecal_QC5 <- metabolomics[, which(metabolomics.sums_fecal > 5)]
# Drop metabolites present in fewer than 3 samples
metabolomics <- dropspc(metabo_fecal_QC5, 3)

# Subset metadata based on sample type (fecal or ileal)
fecal.metadata <- subset(map, Sample_type == "fecal")
ileal.metadata <- subset(map, Sample_type == "ileal")

# Further subset metadata by diet treatment groups
ddgs.fecal.metadata <- subset(fecal.metadata, Diet_trt %in% c("DDGS", "DDGS+AOP"))
rb.fecal.metadata <- subset(fecal.metadata, Diet_trt %in% c("RB", "RB+AOP"))
wm.fecal.metadata <- subset(fecal.metadata, Diet_trt %in% c("WM", "WM+AOP"))

ddgs.ileal.metadata <- subset(ileal.metadata, Diet_trt %in% c("DDGS", "DDGS+AOP"))
rb.ileal.metadata <- subset(ileal.metadata, Diet_trt %in% c("RB", "RB+AOP"))
wm.ileal.metadata <- subset(ileal.metadata, Diet_trt %in% c("WM", "WM+AOP"))

# Extract sample IDs for each group
ddgs.fecal.ids <- rownames(ddgs.fecal.metadata)
rb.fecal.ids <- rownames(rb.fecal.metadata)
wm.fecal.ids <- rownames(wm.fecal.metadata)

ddgs.ileal.ids <- rownames(ddgs.ileal.metadata)
rb.ileal.ids <- rownames(rb.ileal.metadata)
wm.ileal.ids <- rownames(wm.ileal.metadata)

# Subset metabolomics data based on sample IDs
ddgs.fecal.metabo <- metabolomics %>% select(all_of(ddgs.fecal.ids))
rb.fecal.metabo <- metabolomics %>% select(all_of(rb.fecal.ids))
wm.fecal.metabo <- metabolomics %>% select(all_of(wm.fecal.ids))

ddgs.ileal.metabo <- metabolomics %>% select(all_of(ddgs.ileal.ids))
rb.ileal.metabo <- metabolomics %>% select(all_of(rb.ileal.ids))
wm.ileal.metabo <- metabolomics %>% select(all_of(wm.ileal.ids))

# Update metadata to match the selected samples
ddgs.fecal.metadata <- fecal.metadata %>% filter(rownames(fecal.metadata) %in% ddgs.fecal.ids)
rb.fecal.metadata <- fecal.metadata %>% filter(rownames(fecal.metadata) %in% rb.fecal.ids)
wm.fecal.metadata <- fecal.metadata %>% filter(rownames(fecal.metadata) %in% wm.fecal.ids)

ddgs.ileal.metadata <- ileal.metadata %>% filter(rownames(ileal.metadata) %in% ddgs.ileal.ids)
rb.ileal.metadata <- ileal.metadata %>% filter(rownames(ileal.metadata) %in% rb.ileal.ids)
wm.ileal.metadata <- ileal.metadata %>% filter(rownames(ileal.metadata) %in% wm.ileal.ids)

# Define custom colors for plotting
custom_colors <- c(
  "DDGS" = "orange",
  "DDGS+AOP" = "red",
  "RB" = "blue",
  "RB+AOP" = "green",
  "WM" = "#FF66FF",
  "WM+AOP" = "#00FFFF"
)

# ------------------- Maaslin2 and ANCOMBC Analysis -------------------

# DDGS Fecal Metabolomics Analysis
# Prepare metabolic data and metadata
metabolic <- as.data.frame(ddgs.fecal.metabo)
metadata <- ddgs.fecal.metadata

# Create phyloseq object
OTU <- otu_table(as.matrix(metabolic), taxa_are_rows = TRUE)
SAM <- sample_data(metadata)
ps_ddgs_fecal <- phyloseq(OTU, SAM)

# Transform counts to relative abundance
ps.taxa <- transform_sample_counts(ps_ddgs_fecal, function(x) x / sum(x))
ps.taxa.sub <- subset_samples(ps.taxa, Diet_trt %in% c("DDGS", "DDGS+AOP"))

# Run Maaslin2 analysis
mas_1 <- Maaslin2(
  input_data = data.frame(otu_table(ps.taxa.sub)),
  input_metadata = data.frame(sample_data(ps_ddgs_fecal)),
  output = "output_ddgs_fecal",
  min_prevalence = 0,
  fixed_effects = c('Diet_trt'),
  standardize = FALSE
)

# Extract significant results from Maaslin2
mas_res_df <- mas_1$results
fdr_mas <- mas_res_df %>%
  dplyr::filter(qval < 0.05)

# Perform ANCOM-BC analysis
ancombc_results <- ancombc(
  phyloseq = ps.taxa.sub,
  formula = "Diet_trt",
  group = "Diet_trt",
  p_adj_method = "bonferroni",
  struc_zero = TRUE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.01,
  global = TRUE
)

# Extract significant results from ANCOM-BC
res_df <- ancombc_results$res
res <- as.data.frame(res_df)
fdr_ancom <- res %>%
  dplyr::filter(q_val.Diet_trtDDGS.AOP < 0.05)

# Find common significant metabolites between Maaslin2 and ANCOM-BC
common_taxa <- intersect(fdr_mas$feature, fdr_ancom$lfc.taxon)

# Create a dataframe of significant metabolites
significant_taxa <- fdr_mas %>%
  filter(feature %in% common_taxa) %>%
  left_join(fdr_ancom, by = c("feature" = "lfc.taxon"))

# Calculate the score for visualization
significant_taxa$score <- -log10(significant_taxa$qval.x) * sign(significant_taxa$coef.x)

# Assign category based on the sign of the coefficient
significant_taxa$category <- ifelse(significant_taxa$coef.x > 0, "DDGS+AOP", "DDGS")

# Ensure the feature factor levels are ordered by score
significant_taxa$feature <- fct_reorder(significant_taxa$feature, significant_taxa$coef.x)

# Clean up feature names for plotting
significant_taxa$feature <- gsub("\\.", " ", significant_taxa$feature)

# Plot the significant metabolites
ggplot(significant_taxa, aes(x = feature, y = coef.x, fill = category)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  labs(y = "Coefficient", x = "Metabolites") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black"),
    plot.background = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE))

# ------------------- RB Ileal Metabolomics Analysis -------------------

# Prepare metabolic data and metadata
metabolic <- as.data.frame(rb.ileal.metabo)
metadata <- rb.ileal.metadata

# Create phyloseq object
OTU <- otu_table(as.matrix(metabolic), taxa_are_rows = TRUE)
SAM <- sample_data(metadata)
ps_rb_ileal <- phyloseq(OTU, SAM)

# Transform counts to relative abundance
ps.taxa <- transform_sample_counts(ps_rb_ileal, function(x) x / sum(x))
ps.taxa.sub <- subset_samples(ps.taxa, Diet_trt %in% c("RB", "RB+AOP"))

# Run Maaslin2 analysis
mas_1 <- Maaslin2(
  input_data = data.frame(otu_table(ps.taxa.sub)),
  input_metadata = data.frame(sample_data(ps_rb_ileal)),
  output = "output_rb_ileal",
  min_prevalence = 0,
  fixed_effects = c('Diet_trt'),
  standardize = FALSE
)

# Extract significant results from Maaslin2
mas_res_df <- mas_1$results
fdr_mas <- mas_res_df %>%
  dplyr::filter(qval < 0.05)

# Perform ANCOM-BC analysis
ancombc_results <- ancombc(
  phyloseq = ps.taxa.sub,
  formula = "Diet_trt",
  group = "Diet_trt",
  p_adj_method = "bonferroni",
  struc_zero = TRUE,
  neg_lb = TRUE,
  tol = 1e-5,
  max_iter = 100,
  conserve = TRUE,
  alpha = 0.01,
  global = TRUE
)

# Extract significant results from ANCOM-BC
res_df <- ancombc_results$res
res <- as.data.frame(res_df)
fdr_ancom <- res %>%
  dplyr::filter(q_val.Diet_trtRB.AOP < 0.05)

# Find common significant metabolites between Maaslin2 and ANCOM-BC
common_taxa <- intersect(fdr_mas$feature, fdr_ancom$lfc.taxon)

# Create a dataframe of significant metabolites
significant_taxa <- fdr_mas %>%
  filter(feature %in% common_taxa) %>%
  left_join(fdr_ancom, by = c("feature" = "lfc.taxon"))

# Calculate the score for visualization
significant_taxa$score <- -log10(significant_taxa$qval.x) * sign(significant_taxa$coef.x)

# Assign category based on the sign of the coefficient
significant_taxa$category <- ifelse(significant_taxa$coef.x > 0, "RB+AOP", "RB")

# Ensure the feature factor levels are ordered by score
significant_taxa$feature <- fct_reorder(significant_taxa$feature, significant_taxa$coef.x)

# Clean up feature names for plotting
significant_taxa$feature <- gsub("\\.", " ", significant_taxa$feature)

# Plot the significant metabolites
ggplot(significant_taxa, aes(x = feature, y = coef.x, fill = category)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  labs(y = "Coefficient", x = "Metabolites") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = "black"),
    plot.background = element_blank()
  ) +
  guides(fill = guide_legend(reverse = TRUE))

# ------------------- Heatmap of Significant Metabolites -------------------

# Combine significant metabolites from both analyses
ddgs_fecal_sig_metabolic <- as.character(intersect(fdr_mas$feature, fdr_ancom$lfc.taxon))
rb_ileal_sig_metabolic <- as.character(intersect(fdr_mas$feature, fdr_ancom$lfc.taxon))
all_significant_metabolites <- unique(c(ddgs_fecal_sig_metabolic, rb_ileal_sig_metabolic))

# Clean up metabolite names
all_significant_metabolites <- gsub("\\.", " ", all_significant_metabolites)

# Combine phyloseq objects
ps_combined <- merge_phyloseq(ps_ddgs_fecal, ps_rb_ileal)

# Transform counts to relative abundance
ps.taxa.rel <- transform_sample_counts(ps_combined, function(x) x / sum(x) * 100)

# Subset to keep only significant metabolites
ps.taxa.rel.sig <- prune_taxa(all_significant_metabolites, ps.taxa.rel)

# Extract data matrix
matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))

# Update row names
rownames(matrix) <- all_significant_metabolites

# Create annotation data
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))
annotation_col <- data.frame(
  Treatment = metadata_sub$Trt,
  'Sample Type' = metadata_sub$Sample_type,
  Diet = metadata_sub$Diet,
  check.names = FALSE
)
rownames(annotation_col) <- rownames(metadata_sub)

# Define colors for annotations
diet_colors <- c("DDGS" = "#E41A1C", "RB" = "#377EB8", "WM" = "#4DAF4A")
treatment_colors <- c("CTR" = "#984EA3", "AOP" = "#FF7F00")
Sample_Type_colors <- c("fecal" = "#A65628", "ileal" = "#F781BF")

ann_colors <- list(
  Diet = diet_colors,
  Treatment = treatment_colors,
  'Sample Type' = Sample_Type_colors
)

# Scale the rows of the matrix
scaled_matrix <- t(scale(t(matrix)))

# Define color palette for heatmap
color_range <- c(-2, -1, 0, 1, 2)
color_palette <- c("blue", "lightblue", "white", "pink", "red")

# Create the heatmap
Heatmap(
  scaled_matrix,
  name = "Scaled Abundance",
  col = colorRamp2(color_range, color_palette),
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D",
  clustering_distance_columns = "euclidean",
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 8),
  column_title = " ",
  row_title = "Metabolites",
  column_title_gp = gpar(fontsize = 12),
  row_title_gp = gpar(fontsize = 10),
  top_annotation = HeatmapAnnotation(
    df = annotation_col,
    col = ann_colors,
    annotation_name_side = "right",
    annotation_name_gp = gpar(fontsize = 8)
  ),
  column_split = paste(annotation_col$`Sample Type`, annotation_col$Treatment, sep = "_"),
  column_gap = unit(0.4, "mm"),
  heatmap_legend_param = list(
    title = "Scaled Abundance",
    at = color_range,
    labels = c("Low", "", "Medium", "", "High"),
    legend_height = unit(5, "cm"),
    grid_width = unit(0.5, "cm"),
    labels_gp = gpar(fontsize = 8)
  ),
  show_heatmap_legend = TRUE
)
