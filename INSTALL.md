# scQCgenes Installation and Usage Guide

## Installation Steps

### Method 1: Install from Local Directory

```r
# Set working directory to where scQCgenes folder is located
setwd("path/to/directory/containing/scQCgenes")

# Install the package
install.packages("scQCgenes", repos = NULL, type = "source")

# Load the package
library(scQCgenes)
```

### Method 2: Using devtools

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install scQCgenes
devtools::install_local("path/to/scQCgenes")

# Or if you're in the parent directory
devtools::install("scQCgenes")

# Load the package
library(scQCgenes)
```

### Method 3: Using remotes (recommended for deployment)

```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install from local path
remotes::install_local("path/to/scQCgenes", dependencies = TRUE)
```

## Quick Test

After installation, test the package:

```r
library(scQCgenes)

# Test with example gene names
test_genes <- c("mt-Nd1", "mt-Cytb", "Rpl5", "Rps3", 
                "9530026P05Rik", "Gm9968", "Actb", "Gapdh")

# Identify QC genes
qc_genes <- identify_qc_genes(test_genes, species = "mouse")
print(qc_genes)

# Should identify:
# - mt-Nd1, mt-Cytb as mito
# - Rpl5, Rps3 as ribo
# - 9530026P05Rik as rik
# - Gm9968 as gm
```

## Basic Usage Examples

### Example 1: Simple Filtering

```r
library(scQCgenes)
library(Seurat)

# Load your Seurat object
seurat_obj <- readRDS("your_data.rds")

# Filter mitochondrial and ribosomal genes
filtered_obj <- filter_seurat_genes(
  seurat_obj,
  remove_types = c("mito", "ribo")
)

# Check dimensions
cat("Original genes:", nrow(seurat_obj), "\n")
cat("Filtered genes:", nrow(filtered_obj), "\n")
```

### Example 2: Comprehensive QC

```r
# Get detailed information
result <- filter_seurat_genes(
  seurat_obj,
  remove_types = c("mito", "ribo", "rik", "gm"),
  species = "mouse",
  return_removed_genes = TRUE
)

# View summary
print(result$summary)

# View removed genes
head(result$removed_genes)

# Save removed genes list
write.csv(result$removed_genes, "removed_qc_genes.csv", row.names = FALSE)

# Use filtered object
filtered_obj <- result$seurat_obj
```

### Example 3: Visualization

```r
# Before filtering - visualize QC genes
plot_qc_stats(rownames(seurat_obj), species = "mouse")
ggsave("qc_genes_stats.pdf", width = 8, height = 6)

# Get summary statistics
summary_qc_genes(rownames(seurat_obj), species = "mouse")
```

## Integration with Seurat Workflow

```r
library(scQCgenes)
library(Seurat)
library(dplyr)

# Step 1: Load data
seurat_obj <- Read10X(data.dir = "path/to/10x/data")
seurat_obj <- CreateSeuratObject(counts = seurat_obj)

# Step 2: QC visualization (before filtering)
summary_qc_genes(rownames(seurat_obj), species = "mouse")
plot_qc_stats(rownames(seurat_obj), species = "mouse")

# Step 3: Filter QC genes
result <- filter_seurat_genes(
  seurat_obj,
  remove_types = c("mito", "ribo", "rik", "gm"),
  return_removed_genes = TRUE
)

# Step 4: Continue with standard Seurat workflow
seurat_obj <- result$seurat_obj
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Step 5: Find markers (on filtered genes)
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

# Step 6: Select top markers per cluster
library(dplyr)
top5_genes <- markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()
```

## Applying to Your Original Code

Replace your original code:

```r
# OLD CODE:
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
s.genes <- Seurat::cc.genes.updated.2019$s.genes
# ... multiple lines ...
remove_genes <- mito_ribo_genes

# NEW CODE:
library(scQCgenes)
remove_genes <- identify_qc_genes(
  rownames(scedata),
  remove_types = c("mito", "ribo"),
  species = "mouse"
)
```

For filtering FindAllMarkers results:

```r
library(scQCgenes)
library(dplyr)

# After FindAllMarkers
tmpdf <- FindAllMarkers(seurat_obj)
tmpdf$dfpct <- tmpdf$pct.1 - tmpdf$pct.2

# Identify QC genes to remove
qc_genes_to_remove <- identify_qc_genes(
  unique(tmpdf$gene),
  remove_types = c("mito", "ribo", "rik", "gm"),
  species = "mouse"
)

# Filter tmpdf
tmpdf_filtered <- tmpdf %>%
  dplyr::filter(!gene %in% qc_genes_to_remove$Gene)

# Get top 5 per cluster
top5_genes <- tmpdf_filtered %>%
  dplyr::group_by(cluster) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::arrange(desc(dfpct)) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()

# Save results
library(writexl)
write_xlsx(tmpdf_filtered, "FindAllMarkers_filtered.xlsx")
write_xlsx(top5_genes, "Top5_markers_per_cluster.xlsx")
```

## Troubleshooting

### Package Won't Install

```r
# Check dependencies
install.packages(c("stringr", "dplyr", "ggplot2"))

# Try with dependencies
install.packages("scQCgenes", repos = NULL, type = "source", dependencies = TRUE)
```

### Function Not Found

```r
# Make sure package is loaded
library(scQCgenes)

# Check if functions are available
ls("package:scQCgenes")
```

### Seurat Not Available

```r
# Install Seurat
install.packages("Seurat")

# Or from CRAN
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
```

## Getting Help

```r
# View function documentation
?identify_qc_genes
?filter_seurat_genes
?plot_qc_stats

# View vignette
browseVignettes("scQCgenes")
```

## Package Structure

```
scQCgenes/
├── DESCRIPTION              # Package metadata
├── NAMESPACE               # Exported functions
├── LICENSE                 # MIT License
├── README.md              # Package overview
├── INSTALL.md             # This file
├── R/                     # R code
│   ├── identify_qc_genes.R
│   ├── filter_functions.R
│   └── plot_functions.R
└── vignettes/             # Documentation
    └── getting_started.Rmd
```

## Contact

For issues or questions, please contact the package maintainer.
