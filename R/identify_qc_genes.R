#' Identify Genes to Remove for Quality Control
#'
#' Identifies specific types of genes in single-cell RNA-seq data that are
#' typically filtered during quality control, including mitochondrial genes,
#' ribosomal genes, cell cycle genes, and predicted genes (Rik and Gm in mouse).
#'
#' @param gene_names Character vector of gene names (e.g., rownames of Seurat object)
#' @param remove_types Character vector specifying gene types to identify.
#'   Options: "mito" (mitochondrial), "ribo" (ribosomal), "cell_cycle",
#'   "rik" (mouse Rik genes), "gm" (mouse Gm predicted genes).
#'   Default: c("mito", "ribo")
#' @param species Character, either "mouse" or "human". Default: "mouse"
#' @param return_list Logical, if TRUE returns a list with genes grouped by type,
#'   if FALSE returns a data frame. Default: FALSE
#'
#' @return If return_list=FALSE, a data frame with columns:
#'   \itemize{
#'     \item Gene: gene name
#'     \item Type: gene type(s)
#'   }
#'   If return_list=TRUE, a named list where each element contains genes of that type.
#'
#' @export
#' @importFrom stringr str_to_title
#' @examples
#' \dontrun{
#' # Identify mitochondrial and ribosomal genes
#' qc_genes <- identify_qc_genes(rownames(seurat_obj))
#'
#' # Identify all QC gene types
#' qc_genes <- identify_qc_genes(
#'   rownames(seurat_obj),
#'   remove_types = c("mito", "ribo", "cell_cycle", "rik", "gm")
#' )
#'
#' # Return as list
#' qc_list <- identify_qc_genes(rownames(seurat_obj), return_list = TRUE)
#' }
identify_qc_genes <- function(gene_names,
                              remove_types = c("mito", "ribo"),
                              species = "mouse",
                              return_list = FALSE) {
  
  if (!species %in% c("mouse", "human")) {
    stop("Species must be either 'mouse' or 'human'")
  }
  
  valid_types <- c("mito", "ribo", "cell_cycle", "rik", "gm")
  invalid_types <- setdiff(remove_types, valid_types)
  if (length(invalid_types) > 0) {
    stop(paste("Invalid remove_types:", paste(invalid_types, collapse = ", "),
               "\nValid options are:", paste(valid_types, collapse = ", ")))
  }
  
  gene_list <- list()
  
  # Cell cycle genes
  if ("cell_cycle" %in% remove_types) {
    if (requireNamespace("Seurat", quietly = TRUE)) {
      g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
      s.genes <- Seurat::cc.genes.updated.2019$s.genes
      
      if (species == "mouse") {
        g2m.genes <- stringr::str_to_title(g2m.genes)
        s.genes <- stringr::str_to_title(s.genes)
      }
      
      gene_list$cell_cycle <- c(g2m.genes, s.genes)
    } else {
      warning("Seurat package not available, skipping cell cycle genes")
    }
  }
  
  # Mitochondrial genes
  if ("mito" %in% remove_types) {
    if (species == "mouse") {
      mito_pattern <- "^mt-|^Mt-"
    } else {
      mito_pattern <- "^MT-"
    }
    mito_genes <- gene_names[grep(mito_pattern, gene_names)]
    if (length(mito_genes) > 0) {
      gene_list$mito <- mito_genes
    }
  }
  
  # Ribosomal genes
  if ("ribo" %in% remove_types) {
    if (species == "mouse") {
      ribo_pattern <- "^Rp[sl]"
    } else {
      ribo_pattern <- "^RP[SL]"
    }
    ribo_genes <- gene_names[grep(ribo_pattern, gene_names)]
    if (length(ribo_genes) > 0) {
      gene_list$ribo <- ribo_genes
    }
  }
  
  # Rik genes (mouse only)
  if ("rik" %in% remove_types && species == "mouse") {
    rik_genes <- gene_names[grep("Rik[0-9]*$", gene_names)]
    if (length(rik_genes) > 0) {
      gene_list$rik <- rik_genes
    }
  }
  
  # Gm predicted genes (mouse only)
  if ("gm" %in% remove_types && species == "mouse") {
    gm_genes <- gene_names[grep("^Gm[0-9]", gene_names)]
    if (length(gm_genes) > 0) {
      gene_list$gm <- gm_genes
    }
  }
  
  if (return_list) {
    return(gene_list)
  }
  
  # Convert to data frame
  all_genes <- unique(unlist(gene_list))
  
  if (length(all_genes) == 0) {
    return(data.frame(Gene = character(0), Type = character(0), stringsAsFactors = FALSE))
  }
  
  result <- data.frame(
    Gene = all_genes,
    Type = sapply(all_genes, function(g) {
      types <- names(gene_list)[sapply(gene_list, function(x) g %in% x)]
      paste(types, collapse = ", ")
    }),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  return(result)
}
