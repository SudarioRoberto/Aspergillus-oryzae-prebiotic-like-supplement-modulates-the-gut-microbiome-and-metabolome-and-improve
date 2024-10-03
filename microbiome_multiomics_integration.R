# Load necessary libraries for data manipulation, analysis, and visualization
library(caret)
library(tidyverse)
library(BiocManager)
library(phyloseq)
library(Maaslin2)
library(forcats)
library(viridis)
library(knitr)
library(vegan)
library(ape)
library(labdsv)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ANCOMBC)
library(VennDiagram)
library(mixOmics)
library(conflicted)

# Set working directory to the location of your data files
# Please update the path to your specific directory
setwd("/path/to/your/directory")

# Read in the data files: OTU table, taxonomy, and metadata
otu_table <- read.csv("ASV.csv", row.names = 1)
taxonomy <- read.csv("taxon.csv", row.names = 1)
map <- read.csv("metadata.csv", row.names = 1)

# Clean and separate taxonomy information
tax <- taxonomy %>%
  dplyr::select(Taxonomy) %>%
  separate(Taxonomy,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "; ")

# Remove taxonomic prefixes and handle missing data
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

# Replace NA and empty strings with appropriate labels
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean == "__"] <- ""

# Fill in unclassified taxa with higher-level classifications
for (i in 1:nrow(tax.clean)) {
  if (tax.clean[i, 7] != "") {
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i, 2] == "") {
    kingdom <- paste("Unclassified", tax.clean[i, 1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i, 3] == "") {
    phylum <- paste("Unclassified", tax.clean[i, 2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i, 4] == "") {
    class <- paste("Unclassified", tax.clean[i, 3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i, 5] == "") {
    order <- paste("Unclassified", tax.clean[i, 4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i, 6] == "") {
    family <- paste("Unclassified", tax.clean[i, 5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i, 7] == "") {
    tax.clean$Species[i] <- paste("Unclassified", tax.clean$Genus[i], sep = " ")
  }
}

# Overview of the metadata
summary(map)
str(map)

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

# Resolve conflicts in function names
conflict_prefer("select", "dplyr")
conflict_prefer("pca", "mixOmics")
conflicts_prefer(dplyr::filter)

# Subset OTU tables based on sample IDs
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

# Define custom colors for plotting
custom_colors <- c(
  "DDGS" = "orange", 
  "DDGS+AOP" = "red",      
  "RB" = "blue",      
  "RB+AOP" = "green",   
  "WM" = "#FF66FF", 
  "WM+AOP" = "#00FFFF" 
)

# Load metabolomics data
metabolomics <- as.data.frame(t(read.csv("metabolomics.csv", row.names = 1)))

# Subset metabolomics data based on sample IDs
ddgs.fecal.metabo <- metabolomics %>% select(all_of(ddgs.fecal.ids))
rb.fecal.metabo <- metabolomics %>% select(all_of(rb.fecal.ids))
wm.fecal.metabo <- metabolomics %>% select(all_of(wm.fecal.ids))

ddgs.ileal.metabo <- metabolomics %>% select(all_of(ddgs.ileal.ids))
rb.ileal.metabo <- metabolomics %>% select(all_of(rb.ileal.ids))
wm.ileal.metabo <- metabolomics %>% select(all_of(wm.ileal.ids))

# Load and preprocess ATTD data for DDGS fecal samples
ATTD <- read.csv("ddgs_attd.csv", row.names = 1, header = TRUE)
ATTD <- ATTD[, -c(1:5)]  # Remove unwanted columns
ATTD <- ATTD[, 1:6]      # Keep only relevant columns
ATTD <- ATTD[rownames(ATTD) != "", ]  # Remove empty rows

# Create phyloseq object for DDGS fecal samples
OTU <- otu_table(as.matrix(ddgs.fecal.taxa), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(ddgs.fecal.metadata)
ps <- phyloseq(OTU, TAX, SAMPLE)

# Aggregate taxa at the species level
ps.taxa <- tax_glom(ps, taxrank = 'Species', NArm = FALSE)
matrix <- as.data.frame(otu_table(ps.taxa))

# Generate unique row names based on taxonomy
phylum_names <- as.character(tax_table(ps.taxa)[, "Phylum"])
class_names <- as.character(tax_table(ps.taxa)[, "Class"])
order_names <- as.character(tax_table(ps.taxa)[, "Order"])
family_names <- as.character(tax_table(ps.taxa)[, "Family"])
genus_names <- as.character(tax_table(ps.taxa)[, "Genus"])
species_names <- as.character(tax_table(ps.taxa)[, "Species"])

construct_row_name <- function(...) {
  parts <- unlist(list(...))
  non_na_parts <- parts[!is.na(parts) & parts != ""]
  if (length(non_na_parts) < 2) {
    return(non_na_parts[1])
  } else {
    return(paste(tail(non_na_parts, 1), collapse = " "))
  }
}

combined_names <- mapply(construct_row_name, phylum_names, class_names, order_names,
                         family_names, genus_names, species_names)

# Ensure row names are unique
rownames(matrix) <- make.unique(combined_names)
taxa <- as.data.frame(t(matrix))

# Prepare metabolomics data
OTU <- otu_table(as.matrix(ddgs.fecal.metabo), taxa_are_rows = TRUE)
SAMPLE <- sample_data(ddgs.fecal.metadata)
ps <- phyloseq(OTU, SAMPLE)
ps.taxa.rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)
matrix <- data.frame(otu_table(ps.taxa.rel))
metabo <- as.data.frame(t(matrix))
colnames(metabo) <- gsub("\\.", " ", colnames(metabo))
meta <- as.data.frame(sample_data(ps))

# Match ATTD data to metabolomics data
ATTD <- ATTD[match(rownames(metabo), rownames(ATTD)), ]

# Log-transform the data to stabilize variance
metabo <- log(metabo + 0.001)
taxa <- log(taxa + 0.001)
ATTD <- log(ATTD + 0.001)

# Combine data into a list for multi-omics integration
data <- list(metabolic = metabo, taxa = taxa, ATTD = ATTD)
Y <- meta$Diet_trt  # Response variable

# Set seed for reproducibility
set.seed(123)

# Define parameters for PLS models
list.keepX <- c(6, 6)  # Number of variables to keep in X
list.keepY <- c(6, 6)  # Number of variables to keep in Y

# Generate pairwise PLS models
pls1 <- spls(data[["ATTD"]], data[["metabolic"]], keepX = list.keepX, keepY = list.keepY)
pls2 <- spls(data[["metabolic"]], data[["taxa"]], keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(data[["ATTD"]], data[["taxa"]], keepX = list.keepX, keepY = list.keepY)

# Plot variable relationships for each PLS model
plotVar(pls1, cutoff = 0.5, title = "(a) ATTD vs metabolic",
        legend = c("ATTD", "metabolic"), var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2, 2), col = c('darkorchid', 'lightgreen'))

plotVar(pls2, cutoff = 0.5, title = "(b) metabolic vs taxa",
        legend = c("metabolic", "taxa"), var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2, 2), col = c('darkorchid', 'lightgreen'))

plotVar(pls3, cutoff = 0.5, title = "(c) ATTD vs taxa",
        legend = c("ATTD", "taxa"), var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2, 2), col = c('darkorchid', 'lightgreen'))

# Correlation between components
cor(pls1$variates$X, pls1$variates$Y)
cor(pls2$variates$X, pls2$variates$Y)
cor(pls3$variates$X, pls3$variates$Y)

# Create design matrix for DIABLO analysis
design <- matrix(0.1, ncol = length(data), nrow = length(data),
                 dimnames = list(names(data), names(data)))
diag(design) <- 0  # Set diagonal to 0s

# Build initial DIABLO model
basic.diablo.model <- block.splsda(X = data, Y = Y, ncomp = 6, design = design, near.zero.var = TRUE)

# Perform model validation and tuning
perf.diablo <- perf(basic.diablo.model, validation = 'Mfold',
                    folds = 6, nrepeat = 6, near.zero.var = TRUE)
plot(perf.diablo)

# Determine optimal number of components
ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

# Define grid of values for feature selection tuning
test.keepX <- list(
  ATTD = c(5:9, seq(10, 18, 2), seq(20, 30, 5)),
  taxa = c(5:9, seq(10, 18, 2), seq(20, 30, 5)),
  metabolic = c(5:9, seq(10, 18, 2), seq(20, 30, 5))
)

# Run feature selection tuning
tune.TCGA <- tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
                               test.keepX = test.keepX, design = design,
                               validation = 'loo', nrepeat = 6,
                               dist = "centroids.dist", near.zero.var = TRUE)

# Set optimal values of features to retain
list.keepX <- tune.TCGA$choice.keepX

# Build final DIABLO model
final.diablo.model <- block.splsda(X = data, Y = Y, ncomp = ncomp,
                                   keepX = list.keepX, design = design, near.zero.var = TRUE)

# Plot DIABLO results
plotDiablo(final.diablo.model, ncomp = 1, col.per.group = custom_colors[c("DDGS", "DDGS+AOP")])

# Evaluate model performance using AUROC
auroc(final.diablo.model, plot = TRUE, roc.comp = NULL, title = NULL, print = TRUE, roc.block = 3)

# Generate circos plot for variable associations
circosPlot(final.diablo.model,
           cutoff = 0.75,
           line = 1,
           color.blocks = c('darkorchid', 'coral', 'mediumseagreen'),
           color.cor = c("blue", "red"),
           color.Y = custom_colors[c("DDGS", "DDGS+AOP")],
           size.labels = 0.001,
           size.variables = 0.7,
           legend = TRUE,
           legend.title = "",
           linkWidth = 1,
           size.legend = 1,
           showIntraLinks = FALSE,
           var.adj = 0.1,
           block.labels.adj = 0.5)

# Repeat similar steps for RB and WM fecal samples
# ... (Code for RB and WM fecal samples would follow the same structure)

# Repeat similar steps for ileal samples
# ... (Code for ileal samples would follow the same structure)
