# Load required libraries for data analysis and visualization
library(mixOmics)
library(tidyverse)
library(vegan)

# ----------- Data Input -----------

# Read in metabolomics data, taxon data, and metadata
# 'metabolomics.csv' contains metabolomics data with samples as columns and metabolites as rows
# 'taxon.csv' contains taxonomic information for each sample
# 'metadata.csv' contains metadata for each sample, including diet treatment and sample type
metabolomics <- as.data.frame(t(read.csv("metabolomics.csv", row.names = 1)))
taxonomy <- read.csv("taxon.csv", row.names = 1)
map <- read.csv("metadata.csv", row.names = 1)

# ----------- Quality Control for Metabolomics Data -----------

# Calculate the sum of metabolite counts for each sample
metabolomics.sums_fecal <- colSums(metabolomics)

# Filter out metabolites with counts less than or equal to 5
metabo_fecal_QC5 <- metabolomics[, which(metabolomics.sums_fecal > 5)]

# Drop metabolites present in fewer than 3 samples
metabolomics <- dropspc(metabo_fecal_QC5, 3)

# ----------- Subset Metadata by Sample Type -----------

# Subset metadata for fecal and ileal samples
fecal.metadata <- subset(map, Sample_type == "fecal")
ileal.metadata <- subset(map, Sample_type == "ileal")

# Further subset metadata based on different diet treatments for fecal samples
ddgs.fecal.metadata <- subset(fecal.metadata, Diet_trt %in% c("DDGS", "DDGS+AOP"))
rb.fecal.metadata <- subset(fecal.metadata, Diet_trt %in% c("RB", "RB+AOP"))
wm.fecal.metadata <- subset(fecal.metadata, Diet_trt %in% c("WM", "WM+AOP"))

# Subset metadata based on diet treatments for ileal samples
ddgs.ileal.metadata <- subset(ileal.metadata, Diet_trt %in% c("DDGS", "DDGS+AOP"))
rb.ileal.metadata <- subset(ileal.metadata, Diet_trt %in% c("RB", "RB+AOP"))
wm.ileal.metadata <- subset(ileal.metadata, Diet_trt %in% c("WM", "WM+AOP"))

# ----------- Extract Sample IDs -----------

# Extract sample IDs for each group of fecal and ileal samples
ddgs.fecal.ids <- rownames(ddgs.fecal.metadata)
rb.fecal.ids <- rownames(rb.fecal.metadata)
wm.fecal.ids <- rownames(wm.fecal.metadata)

ddgs.ileal.ids <- rownames(ddgs.ileal.metadata)
rb.ileal.ids <- rownames(rb.ileal.metadata)
wm.ileal.ids <- rownames(wm.ileal.metadata)

# Resolve function name conflicts to prioritize dplyr functions
conflicts_prefer(dplyr::filter)

# ----------- Subset Metabolomics Data Based on Sample IDs -----------

# Subset metabolomics data based on fecal sample IDs for each diet group
ddgs.fecal.metabo <- metabolomics %>% select(all_of(ddgs.fecal.ids))
rb.fecal.metabo <- metabolomics %>% select(all_of(rb.fecal.ids))
wm.fecal.metabo <- metabolomics %>% select(all_of(wm.fecal.ids))

# Subset metabolomics data based on ileal sample IDs for each diet group
ddgs.ileal.metabo <- metabolomics %>% select(all_of(ddgs.ileal.ids))
rb.ileal.metabo <- metabolomics %>% select(all_of(rb.ileal.ids))
wm.ileal.metabo <- metabolomics %>% select(all_of(wm.ileal.ids))

# ----------- Update Metadata to Match Selected Samples -----------

# Ensure that metadata corresponds to selected fecal samples
ddgs.fecal.metadata <- fecal.metadata %>% filter(rownames(fecal.metadata) %in% ddgs.fecal.ids)
rb.fecal.metadata <- fecal.metadata %>% filter(rownames(fecal.metadata) %in% rb.fecal.ids)
wm.fecal.metadata <- fecal.metadata %>% filter(rownames(fecal.metadata) %in% wm.fecal.ids)

# Ensure that metadata corresponds to selected ileal samples
ddgs.ileal.metadata <- ileal.metadata %>% filter(rownames(ileal.metadata) %in% ddgs.ileal.ids)
rb.ileal.metadata <- ileal.metadata %>% filter(rownames(ileal.metadata) %in% rb.ileal.ids)
wm.ileal.metadata <- ileal.metadata %>% filter(rownames(ileal.metadata) %in% wm.ileal.ids)

# ----------- Define Custom Colors for Plotting -----------

# Define a custom color palette for the diet groups
custom_colors <- c(
  "DDGS" = "orange",
  "DDGS+AOP" = "red",
  "RB" = "blue",
  "RB+AOP" = "green",
  "WM" = "#FF66FF",
  "WM+AOP" = "#00FFFF"
)

# ----------- DDGS Fecal Metabolomics Analysis -----------

# Prepare metabolomics data and metadata for DDGS fecal samples
metabolic <- as.data.frame(ddgs.fecal.metabo)
metadata <- ddgs.fecal.metadata

# Normalize metabolomics data to relative abundance (i.e., make the sum of each sample equal to 1)
metabolic_normalized <- sweep(metabolic, 2, colSums(metabolic), FUN = "/")

# Log-transform the data (adding a small pseudocount to avoid log(0))
metabolic_log <- log2(metabolic_normalized + 1e-6)

# Transpose the data so that samples are rows and metabolites are columns
metabolic_log_t <- t(metabolic_log)
metabolic_log_t <- as.data.frame(metabolic_log_t)

# Remove near-zero variance variables from the dataset
nzv_indices <- nearZeroVar(metabolic_log_t)

# Ensure that sample IDs in the metadata match the metabolomics data
metadata <- metadata[match(rownames(metabolic_log_t), rownames(metadata)), ]

# Convert the diet treatment variable to a factor for analysis
metadata$Diet_trt <- as.factor(metadata$Diet_trt)

# ----------- Principal Component Analysis (PCA) -----------

# Perform PCA using the mixOmics package
X <- metabolic_log_t  # Data matrix with samples as rows and metabolites as columns

# Set the number of principal components (PCs) to compute
ncomp <- 2

# Remove variables with zero variance from the dataset
zero_var_cols <- which(apply(X, 2, var) == 0)
X_clean <- X[, -zero_var_cols]

# Perform PCA with scaling
pca_model <- pca(X_clean, ncomp = ncomp, scale = TRUE)

# Extract scores (coordinates of samples in PCA space)
scores <- pca_model$variates$X

# Convert PCA scores to a data frame and label components (PC1, PC2, etc.)
scores_df <- as.data.frame(scores)
colnames(scores_df)[1:ncomp] <- paste0("PC", 1:ncomp)
scores_df$Diet_trt <- metadata$Diet_trt  # Add diet treatment as a factor

# ----------- Plotting PCA Results -----------

# Calculate centroids (mean PC values) for each diet treatment group
centroids <- aggregate(scores_df[, c("PC1", "PC2")], 
                       by = list(Diet_trt = scores_df$Diet_trt), 
                       FUN = mean)

# Define axis labels with variance explained by each component
explained_var <- pca_model$prop_expl_var
explained_var_vector <- explained_var$X
var_comp1 <- round(explained_var_vector[1] * 100, 2)  # Variance explained by PC1
var_comp2 <- round(explained_var_vector[2] * 100, 2)  # Variance explained by PC2
xlab_text <- paste0("PC1 (", var_comp1, "%)")
ylab_text <- paste0("PC2 (", var_comp2, "%)")

# Plot the PCA results using ggplot2
p <- ggplot(data = scores_df, aes(x = PC1, y = PC2, color = Diet_trt)) +
  geom_point(size = 4) +  # Plot points for samples
  stat_ellipse(aes(group = Diet_trt), type = "t", level = 0.80) +  # Add 80% confidence ellipses
  scale_color_manual(values = custom_colors) +  # Use custom colors
  geom_segment(data = scores_df, aes(
    xend = centroids[match(Diet_trt, centroids$Diet_trt), "PC1"],
    yend = centroids[match(Diet_trt, centroids$Diet_trt), "PC2"],
    x = PC1, y = PC2), linetype = "dotted") +  # Add dotted lines to centroids
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = c(0.15, 0.89),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7)
  ) +
  labs(x = xlab_text, y = ylab_text) +
  xlim(-12, 12) + ylim(-12, 12)  # Set axis limits

print(p)

# ----------- PERMANOVA on Bray-Curtis Distance -----------

# Normalize data to relative abundance using the total method
metabolic_rel_abund <- decostand(metabolic_log_t, method = "total")

# Calculate Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(metabolic_rel_abund, method = "bray")

# Perform PERMANOVA to test the effect of diet treatment on metabolite profiles
permanova_result <- adonis2(bray_curtis_dist ~ Diet_trt, data = scores_df, permutations = 8000)

# Print PERMANOVA results
print(permanova_result)

print(test.adonis)

# Repeat similar steps for RB and WM ileal metabolomics analyses
# ...

