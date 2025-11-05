#' Filter Seurat Object by Removing QC Genes
#'
#' Removes specified quality control genes from a Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param remove_types Character vector specifying gene types to remove.
#'   Options: "mito", "ribo", "cell_cycle", "rik", "gm".
#'   Default: c("mito", "ribo")
#' @param species Character, either "mouse" or "human". Default: "mouse"
#' @param return_removed_genes Logical, if TRUE returns a list with both
#'   filtered object and removed genes info. Default: FALSE
#'
#' @return If return_removed_genes=FALSE, returns filtered Seurat object.
#'   If return_removed_genes=TRUE, returns a list with:
#'   \itemize{
#'     \item seurat_obj: filtered Seurat object
#'     \item removed_genes: data frame of removed genes
#'     \item n_removed: number of genes removed
#'     \item n_remaining: number of genes remaining
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' # Filter mitochondrial and ribosomal genes
#' filtered_obj <- filter_seurat_genes(seurat_obj)
#'
#' # Filter all QC genes and get info about removed genes
#' result <- filter_seurat_genes(
#'   seurat_obj,
#'   remove_types = c("mito", "ribo", "rik", "gm"),
#'   return_removed_genes = TRUE
#' )
#' filtered_obj <- result$seurat_obj
#' print(paste("Removed", result$n_removed, "genes"))
#' }
filter_seurat_genes <- function(seurat_obj,
                                remove_types = c("mito", "ribo"),
                                species = "mouse",
                                return_removed_genes = FALSE) {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  # Get all gene names
  all_genes <- rownames(seurat_obj)
  
  # Identify genes to remove
  qc_genes <- identify_qc_genes(
    gene_names = all_genes,
    remove_types = remove_types,
    species = species,
    return_list = FALSE
  )
  
  # Genes to keep
  genes_to_keep <- setdiff(all_genes, qc_genes$Gene)
  
  # Filter Seurat object
  filtered_obj <- subset(seurat_obj, features = genes_to_keep)
  
  if (return_removed_genes) {
    return(list(
      seurat_obj = filtered_obj,
      removed_genes = qc_genes,
      n_removed = nrow(qc_genes),
      n_remaining = length(genes_to_keep),
      summary = data.frame(
        Original_genes = length(all_genes),
        Removed_genes = nrow(qc_genes),
        Remaining_genes = length(genes_to_keep),
        Percent_removed = round(nrow(qc_genes) / length(all_genes) * 100, 2)
      )
    ))
  }
  
  return(filtered_obj)
}


#' Filter Gene Expression Matrix by Removing QC Genes
#'
#' Removes specified quality control genes from a gene expression matrix.
#'
#' @param expr_matrix A gene expression matrix with genes as rows
#' @param remove_types Character vector specifying gene types to remove.
#'   Options: "mito", "ribo", "cell_cycle", "rik", "gm".
#'   Default: c("mito", "ribo")
#' @param species Character, either "mouse" or "human". Default: "mouse"
#' @param return_removed_genes Logical, if TRUE returns a list with both
#'   filtered matrix and removed genes info. Default: FALSE
#'
#' @return If return_removed_genes=FALSE, returns filtered matrix.
#'   If return_removed_genes=TRUE, returns a list with filtered matrix and info.
#'
#' @export
#' @examples
#' \dontrun{
#' # Filter expression matrix
#' filtered_matrix <- filter_expression_matrix(expr_matrix)
#' }
filter_expression_matrix <- function(expr_matrix,
                                    remove_types = c("mito", "ribo"),
                                    species = "mouse",
                                    return_removed_genes = FALSE) {
  
  # Get all gene names
  all_genes <- rownames(expr_matrix)
  
  # Identify genes to remove
  qc_genes <- identify_qc_genes(
    gene_names = all_genes,
    remove_types = remove_types,
    species = species,
    return_list = FALSE
  )
  
  # Genes to keep
  genes_to_keep <- setdiff(all_genes, qc_genes$Gene)
  
  # Filter matrix
  filtered_matrix <- expr_matrix[genes_to_keep, , drop = FALSE]
  
  if (return_removed_genes) {
    return(list(
      expr_matrix = filtered_matrix,
      removed_genes = qc_genes,
      n_removed = nrow(qc_genes),
      n_remaining = length(genes_to_keep),
      summary = data.frame(
        Original_genes = length(all_genes),
        Removed_genes = nrow(qc_genes),
        Remaining_genes = length(genes_to_keep),
        Percent_removed = round(nrow(qc_genes) / length(all_genes) * 100, 2)
      )
    ))
  }
  
  return(filtered_matrix)
}
