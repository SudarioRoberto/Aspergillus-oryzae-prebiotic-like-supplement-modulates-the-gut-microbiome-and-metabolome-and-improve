# Load necessary libraries for data manipulation, analysis, and visualization
library(tidyverse)
library(BiocManager)
library(phyloseq)
library(ape)
library(Maaslin2)
library(forcats)
library(viridis)
library(knitr)
library(vegan)
library(labdsv)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ANCOMBC)
library(VennDiagram)

# Set working directory to the location of your data files
# Please update the path to your specific directory
setwd("/path/to/your/directory")

# Read in the data files: OTU table, taxonomy, and metadata
otu_table <- read.csv("ASV.csv", row.names = 1)
taxonomy <- read.csv("taxon.csv", row.names = 1)
map <- read.csv("metadata.csv", row.names = 1)

# ------------------- Quality Control -------------------

# Calculate column sums for the OTU table
asv.sums_fecal <- colSums(otu_table)

# Filter OTUs with counts greater than 5
ASV_fecal_QC5 <- otu_table[, which(asv.sums_fecal > 5)]

# Drop OTUs present in fewer than 3 samples
otu_table <- dropspc(ASV_fecal_QC5, 3)

# ------------------- Clean Taxonomy Names -------------------

# Separate taxonomy into individual levels
tax <- taxonomy %>%
  dplyr::select(Taxonomy) %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order",
                              "Family", "Genus", "Species"), sep = "; ")

# Remove prefixes and handle missing values
tax.clean <- data.frame(
  row.names = row.names(tax),
  Kingdom = str_replace(tax[,1], "k__", ""),
  Phylum = str_replace(tax[,2], "p__", ""),
  Class = str_replace(tax[,3], "c__", ""),
  Order = str_replace(tax[,4], "o__", ""),
  Family = str_replace(tax[,5], "f__", ""),
  Genus = str_replace(tax[,6], "g__", ""),
  Species = str_replace(tax[,7], "s__", ""),
  stringsAsFactors = FALSE
)

# Replace NA and empty strings
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean == "__"] <- ""

# Fill in missing taxonomic levels
for (i in 1:nrow(tax.clean)) {
  if (tax.clean[i,7] != "") {
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == "") {
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == "") {
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == "") {
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == "") {
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == "") {
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == "") {
    tax.clean$Species[i] <- paste("Unclassified", tax.clean$Genus[i], sep = " ")
  }
}

# ------------------- Subset the Data -------------------

# Subset metadata based on sample type
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

# Subset OTU table based on sample IDs
ddgs.fecal.taxa <- otu_table %>% select(all_of(ddgs.fecal.ids))
rb.fecal.taxa <- otu_table %>% select(all_of(rb.fecal.ids))
wm.fecal.taxa <- otu_table %>% select(all_of(wm.fecal.ids))

ddgs.ileal.taxa <- otu_table %>% select(all_of(ddgs.ileal.ids))
rb.ileal.taxa <- otu_table %>% select(all_of(rb.ileal.ids))
wm.ileal.taxa <- otu_table %>% select(all_of(wm.ileal.ids))

# Update metadata to match the selected samples
ddgs.fecal.metadata <- fecal.metadata %>% filter(rownames(fecal.metadata) %in% ddgs.fecal.ids)
rb.fecal.metadata <- fecal.metadata %>% filter(rownames(fecal.metadata) %in% rb.fecal.ids)
wm.fecal.metadata <- fecal.metadata %>% filter(rownames(fecal.metadata) %in% wm.fecal.ids)

ddgs.ileal.metadata <- ileal.metadata %>% filter(rownames(ileal.metadata) %in% ddgs.ileal.ids)
rb.ileal.metadata <- ileal.metadata %>% filter(rownames(ileal.metadata) %in% rb.ileal.ids)
wm.ileal.metadata <- ileal.metadata %>% filter(rownames(ileal.metadata) %in% wm.ileal.ids)

# ------------------- Define Custom Colors -------------------

custom_colors <- c(
  "DDGS" = "orange",
  "DDGS+AOP" = "red",
  "RB" = "blue",
  "RB+AOP" = "green",
  "WM" = "#FF66FF",
  "WM+AOP" = "#00FFFF"
)

# ------------------- Differential Abundance Analysis -------------------

# Function to perform Maaslin2 and ANCOM-BC analyses
perform_differential_analysis <- function(otu_data, tax_data, meta_data, diet_trt_values, diet_trt_label) {
  # Create phyloseq object
  OTU <- otu_table(as.matrix(otu_data), taxa_are_rows = TRUE)
  TAX <- tax_table(as.matrix(tax_data))
  SAMPLE <- sample_data(meta_data)
  ps_object <- phyloseq(OTU, TAX, SAMPLE)
  
  # Transform counts to relative abundance
  ps.taxa.rel <- transform_sample_counts(ps_object, function(x) x / sum(x) * 100)
  
  # Aggregate taxa at the species level
  ps.taxa <- tax_glom(ps_object, taxrank = 'Species', NArm = FALSE)
  ps.taxa.sub <- subset_samples(ps.taxa, Diet_trt %in% diet_trt_values)
  
  # Run Maaslin2 analysis
  mas_res <- Maaslin2(
    input_data = data.frame(otu_table(ps.taxa.sub)),
    input_metadata = data.frame(sample_data(ps_object)),
    output = paste0("output_", diet_trt_label),
    min_prevalence = 0,
    fixed_effects = c('Diet_trt'),
    standardize = FALSE
  )
  
  mas_res_df <- mas_res$results
  fdr_mas <- mas_res_df %>%
    dplyr::filter(qval < 0.05)
  
  # Perform ANCOM-BC analysis
  ancombc_res <- ancombc(
    phyloseq = ps.taxa.sub,
    formula = "Diet_trt",
    lib_cut = 1000,
    group = "Diet_trt",
    p_adj_method = "bonferroni",
    struc_zero = TRUE,
    neg_lb = TRUE,
    tol = 1e-5,
    max_iter = 100,
    conserve = TRUE,
    alpha = 0.05,
    global = TRUE
  )
  
  res_df <- ancombc_res$res
  res <- as.data.frame(res_df)
  
  # Adjust column name for q-values based on diet treatment
  q_val_col <- paste0("q_val.Diet_trt", diet_trt_values[2])
  fdr_ancom <- res %>%
    dplyr::filter(!!sym(q_val_col) < 0.05)
  
  # Find common significant taxa
  common_taxa <- intersect(fdr_mas$feature, fdr_ancom$lfc.taxon)
  
  return(common_taxa)
}

# Perform analysis for each group
ddgs_fecal_sig_taxa <- perform_differential_analysis(
  otu_data = ddgs.fecal.taxa,
  tax_data = tax.clean,
  meta_data = ddgs.fecal.metadata,
  diet_trt_values = c("DDGS", "DDGS+AOP"),
  diet_trt_label = "ddgs_fecal"
)

rb_fecal_sig_taxa <- perform_differential_analysis(
  otu_data = rb.fecal.taxa,
  tax_data = tax.clean,
  meta_data = rb.fecal.metadata,
  diet_trt_values = c("RB", "RB+AOP"),
  diet_trt_label = "rb_fecal"
)

wm_fecal_sig_taxa <- perform_differential_analysis(
  otu_data = wm.fecal.taxa,
  tax_data = tax.clean,
  meta_data = wm.fecal.metadata,
  diet_trt_values = c("WM", "WM+AOP"),
  diet_trt_label = "wm_fecal"
)

ddgs_ileal_sig_taxa <- perform_differential_analysis(
  otu_data = ddgs.ileal.taxa,
  tax_data = tax.clean,
  meta_data = ddgs.ileal.metadata,
  diet_trt_values = c("DDGS", "DDGS+AOP"),
  diet_trt_label = "ddgs_ileal"
)

rb_ileal_sig_taxa <- perform_differential_analysis(
  otu_data = rb.ileal.taxa,
  tax_data = tax.clean,
  meta_data = rb.ileal.metadata,
  diet_trt_values = c("RB", "RB+AOP"),
  diet_trt_label = "rb_ileal"
)

wm_ileal_sig_taxa <- perform_differential_analysis(
  otu_data = wm.ileal.taxa,
  tax_data = tax.clean,
  meta_data = wm.ileal.metadata,
  diet_trt_values = c("WM", "WM+AOP"),
  diet_trt_label = "wm_ileal"
)

# ------------------- Combine Significant Taxa -------------------

# Combine significant taxa from all analyses
all_significant_taxa <- unique(c(
  ddgs_fecal_sig_taxa, rb_fecal_sig_taxa, wm_fecal_sig_taxa,
  ddgs_ileal_sig_taxa, rb_ileal_sig_taxa, wm_ileal_sig_taxa
))

# Remove specific taxa if necessary
all_significant_taxa <- subset(all_significant_taxa, !((all_significant_taxa) %in% c("ASV1140")))

# Merge phyloseq objects
ps_combined <- merge_phyloseq(
  ps_ddgs_fecal, ps_rb_fecal, ps_wm_fecal,
  ps_ddgs_ileal, ps_rb_ileal, ps_wm_ileal
)

# Subset the combined phyloseq object to keep only significant taxa
ps.taxa.rel <- transform_sample_counts(ps_combined, function(x) x / sum(x) * 100)
ps.taxa.rel.sig <- prune_taxa(all_significant_taxa, ps.taxa.rel)

# ------------------- Prepare Data for Heatmap -------------------

# Convert to matrix
matrix <- as.matrix(data.frame(otu_table(ps.taxa.rel.sig)))

# Extract taxonomic levels and create row names
tax_table <- tax_table(ps.taxa.rel.sig)
combined_names <- apply(tax_table, 1, function(x) {
  non_na <- x[!is.na(x) & x != ""]
  if (length(non_na) < 2) return(non_na[1])
  else return(paste(tail(non_na, 1), collapse = " "))
})

rownames(matrix) <- make.unique(combined_names)

# Create annotation data for samples
metadata_sub <- data.frame(sample_data(ps.taxa.rel.sig))
metadata_sub <- metadata_sub %>%
  mutate(
    Trt = ifelse(Trt == "CTR", "- AOP", ifelse(Trt == "AOP", "+ AOP", Trt)),
    Sample_type = ifelse(Sample_type == "fecal", "Fecal", "Ileal")
  )

annotation_col <- data.frame(
  AOP = metadata_sub$Trt,
  'Sample Type' = metadata_sub$Sample_type,
  Diet = metadata_sub$Diet,
  check.names = FALSE
)
rownames(annotation_col) <- rownames(metadata_sub)

# Create annotation data for taxa
annotation_row <- as.factor(tax_table(ps.taxa.rel.sig)[, "Phylum"])
names(annotation_row) <- rownames(matrix)

# Define colors for annotations
diet_colors <- c("DDGS" = "#E41A1C", "RB" = "#377EB8", "WM" = "#4DAF4A")
treatment_colors <- c("- AOP" = "#984EA3", "+ AOP" = "#FF7F00")
Sample_Type_colors <- c("Fecal" = "#A65628", "Ileal" = "#F781BF")
phylum_colors <- colorRampPalette(brewer.pal(8, "Set2"))(length(levels(annotation_row)))
names(phylum_colors) <- levels(annotation_row)

ann_colors <- list(
  Diet = diet_colors,
  AOP = treatment_colors,
  'Sample Type' = Sample_Type_colors,
  Phylum = phylum_colors
)

# Scale rows of the matrix for heatmap
scaled_matrix <- t(scale(t(matrix)))

# ------------------- Create Heatmap -------------------

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
  row_title = "Taxa",
  column_title_gp = gpar(fontsize = 12),
  row_title_gp = gpar(fontsize = 10),
  top_annotation = HeatmapAnnotation(
    df = annotation_col,
    col = ann_colors,
    annotation_name_side = "right",
    annotation_name_gp = gpar(fontsize = 8)
  ),
  left_annotation = rowAnnotation(
    Phylum = annotation_row,
    col = list(Phylum = phylum_colors),
    annotation_name_gp = gpar(fontsize = 8)
  ),
  column_split = paste(annotation_col$AOP, annotation_col$`Sample Type`, sep = "_"),
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
