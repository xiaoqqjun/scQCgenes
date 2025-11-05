# Test Script for scQCgenes Package
# This script tests all major functions of the scQCgenes package

library(scQCgenes)

# Create test data
cat("Creating test gene names...\n")
test_genes <- c(
  # Mitochondrial
  "mt-Nd1", "mt-Nd2", "mt-Cytb", "mt-Co1",
  # Ribosomal
  "Rpl5", "Rpl7", "Rps3", "Rps27",
  # Rik genes
  "9530026P05Rik", "2900026A02Rik", "2310002F09Rik1",
  # Gm genes
  "Gm9968", "Gm17135", "Gm42788",
  # Normal genes
  "Actb", "Gapdh", "Vim", "Krt18",
  # Keep these (lowercase rik or Gm without numbers)
  "Grik2", "Grik3", "Gmap"
)

cat("Total test genes:", length(test_genes), "\n\n")

# Test 1: Identify QC genes
cat("=== Test 1: identify_qc_genes() ===\n")
qc_genes <- identify_qc_genes(
  gene_names = test_genes,
  remove_types = c("mito", "ribo", "rik", "gm"),
  species = "mouse"
)
print(qc_genes)
cat("\nTotal QC genes identified:", nrow(qc_genes), "\n\n")

# Test 2: Return as list
cat("=== Test 2: identify_qc_genes() with return_list=TRUE ===\n")
qc_list <- identify_qc_genes(
  gene_names = test_genes,
  remove_types = c("mito", "ribo", "rik", "gm"),
  species = "mouse",
  return_list = TRUE
)
cat("Mitochondrial genes:", length(qc_list$mito), "\n")
cat("Ribosomal genes:", length(qc_list$ribo), "\n")
cat("Rik genes:", length(qc_list$rik), "\n")
cat("Gm genes:", length(qc_list$gm), "\n\n")

# Test 3: Summary
cat("=== Test 3: summary_qc_genes() ===\n")
summary_qc_genes(
  gene_names = test_genes,
  remove_types = c("mito", "ribo", "rik", "gm"),
  species = "mouse"
)
cat("\n")

# Test 4: Plot (if ggplot2 available)
cat("=== Test 4: plot_qc_stats() ===\n")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p <- plot_qc_stats(
    gene_names = test_genes,
    remove_types = c("mito", "ribo", "rik", "gm"),
    species = "mouse"
  )
  print(p)
  cat("Plot created successfully\n\n")
} else {
  cat("ggplot2 not available, skipping plot test\n\n")
}

# Test 5: Filter expression matrix
cat("=== Test 5: filter_expression_matrix() ===\n")
# Create a simple matrix
test_matrix <- matrix(
  rnorm(length(test_genes) * 10),
  nrow = length(test_genes),
  dimnames = list(test_genes, paste0("Cell", 1:10))
)
cat("Original matrix dimensions:", dim(test_matrix), "\n")

filtered_matrix <- filter_expression_matrix(
  expr_matrix = test_matrix,
  remove_types = c("mito", "ribo", "rik", "gm"),
  species = "mouse"
)
cat("Filtered matrix dimensions:", dim(filtered_matrix), "\n")
cat("Genes removed:", nrow(test_matrix) - nrow(filtered_matrix), "\n\n")

# Test 6: Filter with detailed output
cat("=== Test 6: filter_expression_matrix() with return_removed_genes=TRUE ===\n")
result <- filter_expression_matrix(
  expr_matrix = test_matrix,
  remove_types = c("mito", "ribo", "rik", "gm"),
  species = "mouse",
  return_removed_genes = TRUE
)
print(result$summary)
cat("\n")

# Test 7: Verify correct genes are kept
cat("=== Test 7: Verify filtering logic ===\n")
remaining_genes <- rownames(filtered_matrix)
cat("Remaining genes:", paste(remaining_genes, collapse = ", "), "\n")
cat("\nExpected to remain: Actb, Gapdh, Vim, Krt18, Grik2, Grik3, Gmap\n")
cat("These should be filtered: mt-*, Rpl*, Rps*, *Rik*, Gm[0-9]*\n\n")

# Verify
expected_remaining <- c("Actb", "Gapdh", "Vim", "Krt18", "Grik2", "Grik3", "Gmap")
if (all(expected_remaining %in% remaining_genes) && 
    length(remaining_genes) == length(expected_remaining)) {
  cat("✓ Filtering logic is CORRECT!\n\n")
} else {
  cat("✗ Filtering logic has issues!\n")
  cat("Expected:", paste(expected_remaining, collapse = ", "), "\n")
  cat("Got:", paste(remaining_genes, collapse = ", "), "\n\n")
}

# Test 8: Test different species
cat("=== Test 8: Testing human species ===\n")
human_genes <- c("MT-ND1", "MT-CO1", "RPL5", "RPS3", "ACTB", "GAPDH")
human_qc <- identify_qc_genes(
  gene_names = human_genes,
  remove_types = c("mito", "ribo"),
  species = "human"
)
print(human_qc)
cat("\n")

# Test 9: Cell cycle genes
cat("=== Test 9: Testing cell cycle genes ===\n")
if (requireNamespace("Seurat", quietly = TRUE)) {
  # Get some cell cycle genes (mouse format)
  mouse_cc_genes <- c("Top2a", "Mki67", "Pcna", "Actb")
  cc_qc <- identify_qc_genes(
    gene_names = mouse_cc_genes,
    remove_types = "cell_cycle",
    species = "mouse"
  )
  cat("Cell cycle genes identified:\n")
  print(cc_qc)
} else {
  cat("Seurat not available, skipping cell cycle test\n")
}

cat("\n=== All Tests Completed ===\n")
