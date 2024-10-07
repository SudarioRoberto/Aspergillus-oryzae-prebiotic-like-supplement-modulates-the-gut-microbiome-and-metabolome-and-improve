# Load necessary libraries
library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)


# Load data: OTU table, taxonomy, and metadata
otu_table <- read.csv("ASV.csv", row.names = 1)
taxonomy <- read.csv("taxonomy.csv", row.names = 1)
map <- read.csv("metadata.csv", row.names = 1)

# ----------- Quality Control on OTU Table -----------
# Filter OTUs with counts greater than 5
asv.sums_fecal <- colSums(otu_table)
ASV_fecal_QC5 <- otu_table[, which(asv.sums_fecal > 5)]

# Drop OTUs present in fewer than 3 samples
otu_table <- dropspc(ASV_fecal_QC5, 3)

# ----------- Clean Taxonomy Table -----------
# Separate taxonomy into different levels and clean up prefixes
tax <- taxonomy %>%
  dplyr::select(Taxonomy) %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ")

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__", ""),
                        Phylum = str_replace(tax[,2], "p__", ""),
                        Class = str_replace(tax[,3], "c__", ""),
                        Order = str_replace(tax[,4], "o__", ""),
                        Family = str_replace(tax[,5], "f__", ""),
                        Genus = str_replace(tax[,6], "g__", ""),
                        Species = str_replace(tax[,7], "s__", ""),
                        stringsAsFactors = FALSE)

# Handle missing or unclassified taxa
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean == "__"] <- ""

for (i in 1:nrow(tax.clean)) {
  if (tax.clean[i,7] != "") {
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else {
    tax.clean[i, 2:7] <- apply(tax.clean[i, 2:7], 2, function(x) ifelse(x == "", paste0("Unclassified ", tax.clean[i, which.max(x != "")]), x))
  }
}

# ----------- Subset Metadata Based on Sample Type -----------
# Subset metadata for fecal and ileal samples
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

# ----------- Subset OTU Table Based on Sample IDs -----------
# Subset OTU table for each diet group
ddgs.fecal.taxa <- otu_table %>% select(all_of(ddgs.fecal.ids))
rb.fecal.taxa <- otu_table %>% select(all_of(rb.fecal.ids))
wm.fecal.taxa <- otu_table %>% select(all_of(wm.fecal.ids))

ddgs.ileal.taxa <- otu_table %>% select(all_of(ddgs.ileal.ids))
rb.ileal.taxa <- otu_table %>% select(all_of(rb.ileal.ids))
wm.ileal.taxa <- otu_table %>% select(all_of(wm.ileal.ids))

# ----------- Custom Colors for Plotting -----------
custom_colors <- c(
  "DDGS" = "orange", 
  "DDGS+AOP" = "red",      
  "RB" = "blue",      
  "RB+AOP" = "green",   
  "WM" = "#FF66FF", 
  "WM+AOP" = "#00FFFF" 
)

# ----------- Alpha Diversity Analysis -----------
# Create Phyloseq object for DDGS fecal samples
OTU <- otu_table(as.matrix(ddgs.fecal.taxa), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(ddgs.fecal.metadata)
ps <- phyloseq(OTU, TAX, SAMPLE)

# Plot richness (Shannon index)
plot_richness(ps, x = "Diet_trt", measures = c("Shannon")) +
  geom_boxplot(aes(fill = Diet_trt), outlier.shape = NA) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  labs(x = "Diet", y = "Shannon's H", title = NULL) +
  coord_cartesian(ylim = c(4, 5.5))

# Estimate richness and perform pairwise Wilcoxon test
rich <- estimate_richness(ps, measures = c("Shannon"))
wilcox.observed <- pairwise.wilcox.test(rich$Shannon, sample_data(ps)$Diet_trt)

# ----------- Beta Diversity Analysis -----------
# Normalize OTU counts and perform ordination (PCoA with Bray-Curtis distance)
ps <- transform_sample_counts(ps, function(x) x / sum(x))
ordination <- ordinate(ps, method = "PCoA", distance = "bray")

# Plot PCoA results
ordination_scores <- ordination$vectors[, 1:2]
explained_variance <- ordination$values$Relative_eig * 100
xlab_text <- paste0("PCo 1 (", round(explained_variance[1], 2), "%)")
ylab_text <- paste0("PCo 2 (", round(explained_variance[2], 2), "%)")

ordination_df <- as.data.frame(ordination_scores)
ordination_df$SampleID <- rownames(ordination_df)
ordination_df$Diet_trt <- sample_data(ps)$Diet_trt[match(ordination_df$SampleID, rownames(sample_data(ps)))]

centroids <- aggregate(ordination_df[, 1:2], list(Diet_trt = ordination_df$Diet_trt), mean)
names(centroids)[2:3] <- c("PCo 1", "PCo 2")

# PCoA plot
p <- ggplot(data = ordination_df, aes(x = Axis.1, y = Axis.2, color = Diet_trt)) +
  geom_point(size = 4) +
  stat_ellipse(aes(group = Diet_trt), type = "t", level = 0.45) +
  scale_color_manual(values = custom_colors) +
  theme_classic() +
  labs(x = xlab_text, y = ylab_text) +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5))

print(p)

# ----------- PERMANOVA Analysis -----------
# Perform PERMANOVA on Bray-Curtis distance matrix
metadata <- data.frame(sample_data(ps))
dist <- distance(ps, method = "bray")
test.adonis <- adonis2(dist ~ Diet_trt, data = metadata)

# Print PERMANOVA result
print(test.adonis)
