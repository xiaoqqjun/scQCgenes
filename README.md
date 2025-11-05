# scQCgenes: Quality Control Gene Filtering for Single-Cell RNA-seq

`scQCgenes` is an R package that identifies and filters quality control genes in single-cell RNA-seq data. It helps remove genes that may confound downstream analysis, including:

- **Mitochondrial genes** (mt- or MT-)
- **Ribosomal genes** (Rp[sl] or RP[SL])
- **Cell cycle genes** (from Seurat's updated list)
- **Predicted genes** (mouse-specific):
  - Rik genes (RIKEN cDNA clones)
  - Gm genes (predicted genes with Gm + numbers)

## Installation

```r
devtools::install_github("xiaoqqjun/scQCgenes")
```

## Quick Start

### 1. Identify QC Genes

```r
library(scQCgenes)

# Load your Seurat object
library(Seurat)
# seurat_obj <- readRDS("your_seurat_object.rds")

# Identify mitochondrial and ribosomal genes
qc_genes <- identify_qc_genes(
  rownames(seurat_obj),
  remove_types = c("mito", "ribo")
)

head(qc_genes)
#   Gene    Type
# 1 mt-Nd1  mito
# 2 mt-Nd2  mito
# 3 Rpl5    ribo
# 4 Rps3    ribo
```

### 2. Filter Seurat Object

```r
# Remove mitochondrial and ribosomal genes
filtered_obj <- filter_seurat_genes(
  seurat_obj,
  remove_types = c("mito", "ribo")
)

# Remove all QC gene types (for mouse)
filtered_obj <- filter_seurat_genes(
  seurat_obj,
  remove_types = c("mito", "ribo", "cell_cycle", "rik", "gm"),
  species = "mouse"
)

# Get detailed information about removed genes
result <- filter_seurat_genes(
  seurat_obj,
  remove_types = c("mito", "ribo", "rik", "gm"),
  return_removed_genes = TRUE
)

# Access filtered object
filtered_obj <- result$seurat_obj

# Check summary
print(result$summary)
#   Original_genes Removed_genes Remaining_genes Percent_removed
# 1          20000          1500           18500            7.5
```

### 3. Visualize QC Genes

```r
# Plot gene counts by type
plot_qc_stats(
  rownames(seurat_obj),
  remove_types = c("mito", "ribo", "cell_cycle", "rik", "gm")
)

# Plot as percentages
plot_qc_stats(
  rownames(seurat_obj),
  plot_type = "percentage"
)

# Get summary table
summary_qc_genes(
  rownames(seurat_obj),
  remove_types = c("mito", "ribo", "cell_cycle", "rik", "gm")
)
```

### 4. Filter Expression Matrix

```r
# If you have a raw expression matrix
filtered_matrix <- filter_expression_matrix(
  expr_matrix,
  remove_types = c("mito", "ribo")
)
```

## Examples

### Mouse Single-Cell Analysis

```r
library(scQCgenes)
library(Seurat)

# Load mouse scRNA-seq data
seurat_obj <- readRDS("mouse_scRNA.rds")

# Remove all QC genes
result <- filter_seurat_genes(
  seurat_obj,
  remove_types = c("mito", "ribo", "cell_cycle", "rik", "gm"),
  species = "mouse",
  return_removed_genes = TRUE
)

# Check what was removed
summary_qc_genes(rownames(seurat_obj), species = "mouse")

# Visualize
plot_qc_stats(rownames(seurat_obj), species = "mouse")

# Continue with downstream analysis
filtered_obj <- result$seurat_obj
filtered_obj <- NormalizeData(filtered_obj)
filtered_obj <- FindVariableFeatures(filtered_obj)
# ... etc
```

### Human Single-Cell Analysis

```r
# For human data
result <- filter_seurat_genes(
  human_seurat_obj,
  remove_types = c("mito", "ribo", "cell_cycle"),
  species = "human",
  return_removed_genes = TRUE
)

# Note: rik and gm types are mouse-specific and ignored for human
```

### Custom Filtering

```r
# Only remove mitochondrial genes
qc_genes <- identify_qc_genes(
  rownames(seurat_obj),
  remove_types = "mito"
)

# Remove multiple types
qc_genes <- identify_qc_genes(
  rownames(seurat_obj),
  remove_types = c("mito", "ribo", "rik", "gm")
)

# Get as a list grouped by type
qc_list <- identify_qc_genes(
  rownames(seurat_obj),
  remove_types = c("mito", "ribo", "rik", "gm"),
  return_list = TRUE
)

# Access specific types
mito_genes <- qc_list$mito
rik_genes <- qc_list$rik
```

## Gene Types

### Mitochondrial Genes
- **Mouse**: Genes starting with `mt-` or `Mt-`
- **Human**: Genes starting with `MT-`
- Examples: mt-Nd1, mt-Cytb, MT-CO1

### Ribosomal Genes
- **Mouse**: Genes starting with `Rpl` or `Rps`
- **Human**: Genes starting with `RPL` or `RPS`
- Examples: Rpl5, Rps3, RPL7, RPS27

### Cell Cycle Genes
- S phase and G2M phase genes from Seurat's updated 2019 list
- Automatically converted to proper case for mouse (Title case)
- Examples: Top2a, Mki67, PCNA

### Rik Genes (Mouse Only)
- RIKEN cDNA clones, predicted genes with unknown function
- Pattern: Ends with `Rik` or `Rik` + digits
- Examples: 9530026P05Rik, 2900026A02Rik, 2310002F09Rik1

### Gm Genes (Mouse Only)
- Predicted genes in mouse genome
- Pattern: Starts with `Gm` followed by digits
- Examples: Gm9968, Gm17135, Gm42788

**Note**: Genes like `Grik2` (lowercase rik) and `Gmap` (Gm not followed by digits) are **NOT** filtered.

## Functions

| Function | Description |
|----------|-------------|
| `identify_qc_genes()` | Identify QC genes from gene names |
| `filter_seurat_genes()` | Filter Seurat object by removing QC genes |
| `filter_expression_matrix()` | Filter expression matrix by removing QC genes |
| `plot_qc_stats()` | Visualize QC gene statistics |
| `summary_qc_genes()` | Print summary table of QC genes |

## Parameters

### `remove_types`
Character vector specifying which gene types to identify/remove:
- `"mito"`: Mitochondrial genes
- `"ribo"`: Ribosomal genes
- `"cell_cycle"`: Cell cycle genes
- `"rik"`: Rik genes (mouse only)
- `"gm"`: Gm predicted genes (mouse only)

### `species`
- `"mouse"`: For mouse data (default)
- `"human"`: For human data

## Why Filter These Genes?

1. **Mitochondrial genes**: High expression may indicate cell stress or poor quality
2. **Ribosomal genes**: Can mask biological variation in clustering
3. **Cell cycle genes**: May confound cell type identification
4. **Rik/Gm genes**: Predicted genes with unknown function, difficult to interpret

## Citation

If you use `scQCgenes` in your research, please cite:

```
scQCgenes: Quality Control Gene Filtering for Single-Cell RNA-seq
Version 0.1.0
https://github.com/xiaoqqjun/scQCgenes
```

## License

MIT License

## Contact

For questions or issues, please open an issue on GitHub.

## See Also

- [Seurat](https://satijalab.org/seurat/)
- [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/)
- [scater](https://bioconductor.org/packages/scater/)
