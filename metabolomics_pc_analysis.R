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
library(ggrepel)

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

# Resolve function name conflicts
conflicts_prefer(dplyr::filter)

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

# --------------- DDGS Fecal Metabolomics Analysis ------------------

# Prepare metabolic data and metadata
metabolic <- as.data.frame(ddgs.fecal.metabo)
metadata <- ddgs.fecal.metadata

# Create phyloseq object
OTU <- otu_table(as.matrix(metabolic), taxa_are_rows = TRUE)
SAM <- sample_data(metadata)
ps <- phyloseq(OTU, SAM)

# Normalize counts to relative abundance
ps <- transform_sample_counts(ps, function(x) x / sum(x))

# Perform ordination using Principal Coordinates Analysis (PCoA) with Bray-Curtis distance
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Extract ordination scores (first two axes)
ordination_scores <- ordination$vectors[, 1:2]

# Calculate percentage of variance explained by each axis
explained_variance <- ordination$values$Relative_eig * 100
xlab_text <- paste0("PC1 (", round(explained_variance[1], 2), "%)")
ylab_text <- paste0("PC2 (", round(explained_variance[2], 2), "%)")

# Prepare data frame for plotting and environmental fitting
ordination_df <- as.data.frame(ordination_scores)
ordination_df$SampleID <- rownames(ordination_df)
ordination_df$Diet_trt <- sample_data(ps)$Diet_trt[match(ordination_df$SampleID, rownames(sample_data(ps)))]

# Transpose and normalize metabolic data for environmental fitting
metabolic_transposed <- as.data.frame(t(metabolic))
metabolic_normalized <- decostand(metabolic_transposed, method = "total") * 100

# Fit environmental variables (metabolites) onto ordination
env_fit <- envfit(ordination_df[, 1:2], metabolic_normalized, perm = 999)

# Extract vectors and p-values
vectors_data <- as.data.frame(env_fit$vectors$arrows)
p_values <- env_fit$vectors$pvals

# Combine vectors data and p-values
vectors_data$p_value <- p_values

# Filter significant vectors with p-value < 0.01
significant_vectors <- vectors_data[vectors_data$p_value < 0.01, ]
significant_vectors$label <- rownames(significant_vectors)
significant_vectors$label <- gsub("\\.", " ", significant_vectors$label)

# Calculate centroids for spider plot visualization
centroids <- aggregate(ordination_df[, 1:2], list(Diet_trt = ordination_df$Diet_trt), mean)
names(centroids)[2:3] <- c("PC1", "PC2")

# Plot ordination with environmental vectors and centroids
p <- ggplot(data = ordination_df, aes(x = Axis.1, y = Axis.2, color = Diet_trt)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Diet_trt), type = "t", level = 0.80) +
  scale_color_manual(values = custom_colors) +
  # Draw lines from samples to centroids
  geom_segment(data = ordination_df, aes(
    xend = centroids[match(Diet_trt, centroids$Diet_trt), "PC1"],
    yend = centroids[match(Diet_trt, centroids$Diet_trt), "PC2"],
    x = Axis.1, y = Axis.2), linetype = "dotted") +
  # Add axes lines
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Add significant environmental vectors
  geom_segment(data = significant_vectors, aes(
    x = 0, y = 0, xend = Axis.1 / 2.8, yend = Axis.2 / 2.8),
    arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  # Label significant vectors
  geom_text_repel(data = significant_vectors, aes(
    x = Axis.1 / 2.3, y = Axis.2 / 2.3, label = label),
    color = "black", size = 2.8) +
  # Customize theme and labels
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = c(0.80, 0.89),
    legend.direction = "vertical",
    legend.box.just = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.key = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.7)
  ) +
  labs(x = xlab_text, y = ylab_text) +
  scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0.02, 0.02))

# Display the plot
print(p)

# Perform PERMANOVA to test for significant differences between diet treatments
metadata_df <- data.frame(sample_data(ps))
dist <- distance(ps, method = "bray")
test.adonis <- adonis2(dist ~ Diet_trt, data = metadata_df)

# Print PERMANOVA results
print(test.adonis)

# --------------- RB Fecal Metabolomics Analysis ------------------

# Prepare metabolic data and metadata
metabolic <- as.data.frame(rb.fecal.metabo)
metadata <- rb.fecal.metadata

# Create phyloseq object
OTU <- otu_table(as.matrix(metabolic), taxa_are_rows = TRUE)
SAM <- sample_data(metadata)
ps <- phyloseq(OTU, SAM)

# Optional: Remove specific samples if needed (commented out)
# ps <- prune_samples(sample_names(ps) != "SampleID_to_remove", ps)

# Normalize counts to relative abundance
ps <- transform_sample_counts(ps, function(x) x / sum(x))

# Perform ordination using PCoA with Bray-Curtis distance
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Extract ordination scores and explained variance
ordination_scores <- ordination$vectors[, 1:2]
explained_variance <- ordination$values$Relative_eig * 100
xlab_text <- paste0("PC1 (", round(explained_variance[1], 2), "%)")
ylab_text <- paste0("PC2 (", round(explained_variance[2], 2), "%)")

# Prepare data frame for plotting
ordination_df <- as.data.frame(ordination_scores)
ordination_df$SampleID <- rownames(ordination_df)
ordination_df$Diet_trt <- sample_data(ps)$Diet_trt[match(ordination_df$SampleID, rownames(sample_data(ps)))]

# Calculate centroids for spider plot visualization
centroids <- aggregate(ordination_df[, 1:2], list(Diet_trt = ordination_df$Diet_trt), mean)
names(centroids)[2:3] <- c("PC1", "PC2")

# Plot ordination
p <- ggplot(data = ordination_df, aes(x = Axis.1, y = Axis.2, color = Diet_trt)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Diet_trt), type = "t", level = 0.45) +
  scale_color_manual(values = custom_colors) +
  geom_segment(data = ordination_df, aes(
    xend = centroids[match(Diet_trt, centroids$Diet_trt), "PC1"],
    yend = centroids[match(Diet_trt, centroids$Diet_trt), "PC2"],
    x = Axis.1, y = Axis.2), linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = c(0.80, 0.89),
    legend.direction = "vertical",
    legend.box.just = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.key = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.7)
  ) +
  labs(x = xlab_text, y = ylab_text) +
  scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0.02, 0.02))

# Display the plot
print(p)

# Perform PERMANOVA
metadata_df <- data.frame(sample_data(ps))
dist <- distance(ps, method = "bray")
test.adonis <- adonis2(dist ~ Diet_trt, data = metadata_df)

# Print PERMANOVA results
print(test.adonis)

# --------------- WM Fecal Metabolomics Analysis ------------------

# Prepare metabolic data and metadata
metabolic <- as.data.frame(wm.fecal.metabo)
metadata <- wm.fecal.metadata

# Create phyloseq object
OTU <- otu_table(as.matrix(metabolic), taxa_are_rows = TRUE)
SAM <- sample_data(metadata)
ps <- phyloseq(OTU, SAM)

# Normalize counts to relative abundance
ps <- transform_sample_counts(ps, function(x) x / sum(x))

# Perform ordination using PCoA with Bray-Curtis distance
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Extract ordination scores and explained variance
ordination_scores <- ordination$vectors[, 1:2]
explained_variance <- ordination$values$Relative_eig * 100
xlab_text <- paste0("PC1 (", round(explained_variance[1], 2), "%)")
ylab_text <- paste0("PC2 (", round(explained_variance[2], 2), "%)")

# Prepare data frame for plotting and environmental fitting
ordination_df <- as.data.frame(ordination_scores)
ordination_df$SampleID <- rownames(ordination_df)
ordination_df$Diet_trt <- sample_data(ps)$Diet_trt[match(ordination_df$SampleID, rownames(sample_data(ps)))]

# Transpose and normalize metabolic data for environmental fitting
metabolic_transposed <- as.data.frame(t(metabolic))
metabolic_normalized <- decostand(metabolic_transposed, method = "total") * 100

# Fit environmental variables onto ordination
env_fit <- envfit(ordination_df[, 1:2], metabolic_normalized, perm = 999)

# Extract significant vectors
vectors_data <- as.data.frame(env_fit$vectors$arrows)
p_values <- env_fit$vectors$pvals
vectors_data$p_value <- p_values
significant_vectors <- vectors_data[vectors_data$p_value < 0.01, ]
significant_vectors$label <- rownames(significant_vectors)
significant_vectors$label <- gsub("\\.", " ", significant_vectors$label)

# Calculate centroids
centroids <- aggregate(ordination_df[, 1:2], list(Diet_trt = ordination_df$Diet_trt), mean)
names(centroids)[2:3] <- c("PC1", "PC2")

# Plot ordination
p <- ggplot(data = ordination_df, aes(x = Axis.1, y = Axis.2, color = Diet_trt)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Diet_trt), type = "t", level = 0.80) +
  scale_color_manual(values = custom_colors) +
  geom_segment(data = ordination_df, aes(
    xend = centroids[match(Diet_trt, centroids$Diet_trt), "PC1"],
    yend = centroids[match(Diet_trt, centroids$Diet_trt), "PC2"],
    x = Axis.1, y = Axis.2), linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_segment(data = significant_vectors, aes(
    x = 0, y = 0, xend = Axis.1 / 2.8, yend = Axis.2 / 2.8),
    arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  geom_text_repel(data = significant_vectors, aes(
    x = Axis.1 / 2.3, y = Axis.2 / 2.3, label = label),
    color = "black", size = 3) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = c(0.82, 0.89),
    legend.direction = "vertical",
    legend.box.just = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.key = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.7)
  ) +
  labs(x = xlab_text, y = ylab_text) +
  scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0.02, 0.02))

# Display the plot
print(p)

# Perform PERMANOVA
metadata_df <- data.frame(sample_data(ps))
dist <- distance(ps, method = "bray")
test.adonis <- adonis2(dist ~ Diet_trt, data = metadata_df)

# Print PERMANOVA results
print(test.adonis)

# --------------- DDGS Ileal Metabolomics Analysis ------------------

# Prepare metabolic data and metadata
metabolic <- as.data.frame(ddgs.ileal.metabo)
metadata <- ddgs.ileal.metadata

# Create phyloseq object
OTU <- otu_table(as.matrix(metabolic), taxa_are_rows = TRUE)
SAM <- sample_data(metadata)
ps <- phyloseq(OTU, SAM)

# Normalize counts to relative abundance
ps <- transform_sample_counts(ps, function(x) x / sum(x))

# Perform ordination using PCoA with Bray-Curtis distance
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Extract ordination scores and explained variance
ordination_scores <- ordination$vectors[, 1:2]
explained_variance <- ordination$values$Relative_eig * 100
xlab_text <- paste0("PC1 (", round(explained_variance[1], 2), "%)")
ylab_text <- paste0("PC2 (", round(explained_variance[2], 2), "%)")

# Prepare data frame for plotting
ordination_df <- as.data.frame(ordination_scores)
ordination_df$SampleID <- rownames(ordination_df)
ordination_df$Diet_trt <- sample_data(ps)$Diet_trt[match(ordination_df$SampleID, rownames(sample_data(ps)))]

# Calculate centroids
centroids <- aggregate(ordination_df[, 1:2], list(Diet_trt = ordination_df$Diet_trt), mean)
names(centroids)[2:3] <- c("PC1", "PC2")

# Plot ordination
p <- ggplot(data = ordination_df, aes(x = Axis.1, y = Axis.2, color = Diet_trt)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Diet_trt), type = "t", level = 0.45) +
  scale_color_manual(values = custom_colors) +
  geom_segment(data = ordination_df, aes(
    xend = centroids[match(Diet_trt, centroids$Diet_trt), "PC1"],
    yend = centroids[match(Diet_trt, centroids$Diet_trt), "PC2"],
    x = Axis.1, y = Axis.2), linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = c(0.80, 0.89),
    legend.direction = "vertical",
    legend.box.just = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.key = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.7)
  ) +
  labs(x = xlab_text, y = ylab_text) +
  scale_x_continuous(limits = c(-0.5, 0.55), expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0.02, 0.02))

# Display the plot
print(p)

# Perform PERMANOVA
metadata_df <- data.frame(sample_data(ps))
dist <- distance(ps, method = "bray")
test.adonis <- adonis2(dist ~ Diet_trt, data = metadata_df)

# Print PERMANOVA results
print(test.adonis)

# Repeat similar steps for RB and WM ileal metabolomics analyses
# ...

