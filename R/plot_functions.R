#' Plot QC Gene Statistics
#'
#' Creates a bar plot showing the number and percentage of each type of QC gene.
#'
#' @param gene_names Character vector of gene names
#' @param remove_types Character vector specifying gene types to plot.
#'   Default: c("mito", "ribo", "cell_cycle", "rik", "gm")
#' @param species Character, either "mouse" or "human". Default: "mouse"
#' @param plot_type Character, either "count" or "percentage". Default: "count"
#'
#' @return A ggplot2 object
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_bar geom_text theme_minimal labs coord_flip
#' @importFrom dplyr group_by summarise n
#' @examples
#' \dontrun{
#' # Plot QC gene statistics
#' plot_qc_stats(rownames(seurat_obj))
#'
#' # Plot as percentages
#' plot_qc_stats(rownames(seurat_obj), plot_type = "percentage")
#' }
plot_qc_stats <- function(gene_names,
                         remove_types = c("mito", "ribo", "cell_cycle", "rik", "gm"),
                         species = "mouse",
                         plot_type = "count") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for plotting")
  }
  
  # Get QC genes by type
  qc_list <- identify_qc_genes(
    gene_names = gene_names,
    remove_types = remove_types,
    species = species,
    return_list = TRUE
  )
  
  # Count genes by type
  stats_df <- data.frame(
    Type = names(qc_list),
    Count = sapply(qc_list, length),
    stringsAsFactors = FALSE
  )
  
  # Calculate percentages
  total_genes <- length(gene_names)
  stats_df$Percentage <- round(stats_df$Count / total_genes * 100, 2)
  
  # Reorder by count
  stats_df <- stats_df[order(stats_df$Count, decreasing = TRUE), ]
  stats_df$Type <- factor(stats_df$Type, levels = stats_df$Type)
  
  # Create plot
  if (plot_type == "percentage") {
    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = Type, y = Percentage)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = paste0(Percentage, "%")),
                        vjust = -0.5, size = 4) +
      ggplot2::labs(
        title = "QC Genes by Type (Percentage)",
        x = "Gene Type",
        y = "Percentage of Total Genes (%)"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  } else {
    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = Type, y = Count)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = Count),
                        vjust = -0.5, size = 4) +
      ggplot2::labs(
        title = "QC Genes by Type (Count)",
        x = "Gene Type",
        y = "Number of Genes"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  }
  
  return(p)
}


#' Summary of QC Genes
#'
#' Prints a summary table of QC genes identified in the dataset.
#'
#' @param gene_names Character vector of gene names
#' @param remove_types Character vector specifying gene types to summarize.
#'   Default: c("mito", "ribo", "cell_cycle", "rik", "gm")
#' @param species Character, either "mouse" or "human". Default: "mouse"
#'
#' @return A data frame with summary statistics
#'
#' @export
#' @examples
#' \dontrun{
#' # Get summary of QC genes
#' summary_qc_genes(rownames(seurat_obj))
#' }
summary_qc_genes <- function(gene_names,
                            remove_types = c("mito", "ribo", "cell_cycle", "rik", "gm"),
                            species = "mouse") {
  
  # Get QC genes by type
  qc_list <- identify_qc_genes(
    gene_names = gene_names,
    remove_types = remove_types,
    species = species,
    return_list = TRUE
  )
  
  # Create summary
  total_genes <- length(gene_names)
  
  summary_df <- data.frame(
    Gene_Type = names(qc_list),
    Count = sapply(qc_list, length),
    Percentage = round(sapply(qc_list, length) / total_genes * 100, 2),
    stringsAsFactors = FALSE
  )
  
  # Add total row
  total_unique <- length(unique(unlist(qc_list)))
  summary_df <- rbind(
    summary_df,
    data.frame(
      Gene_Type = "Total (unique)",
      Count = total_unique,
      Percentage = round(total_unique / total_genes * 100, 2)
    )
  )
  
  # Add remaining genes
  summary_df <- rbind(
    summary_df,
    data.frame(
      Gene_Type = "Remaining genes",
      Count = total_genes - total_unique,
      Percentage = round((total_genes - total_unique) / total_genes * 100, 2)
    )
  )
  
  cat("\n=== QC Genes Summary ===\n")
  cat(paste("Total genes:", total_genes, "\n"))
  cat(paste("Species:", species, "\n\n"))
  
  print(summary_df, row.names = FALSE)
  
  return(invisible(summary_df))
}
