#' Merge results from multiple ActivePathways analyses
#'
#' This function combines results from multiple ActivePathways analyses into a single set of files
#' for visualization in Cytoscape. This is particularly useful for comparing results with and without
#' directional penalties.
#'
#' @param enriched_pathways A data.table returned by ActivePathways
#' @param enriched_pathways_directional A data.table returned by ActivePathways
#' @param gmt_file Path to GMT file
#' @param output_prefix A string prefix for output files
#' @param col_colors A character vector of colors for each test (must match length of tests)
#' @param tests A character vector of names for the data sources (e.g., c('rna', 'protein', 'combined')) or NULL
#'
#' @return A list containing the merged results
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' # Run two different ActivePathways analyses
#' enriched_pathways <- ActivePathways(
#'   pval_matrix, gmt = fname_GMT2, cytoscape_file_tag = "original_")
#'
#' enriched_pathways_directional <- ActivePathways(
#'   pval_matrix, gmt = fname_GMT2, cytoscape_file_tag = "directional_",
#'   merge_method = "DPM", scores_direction = dir_matrix, 
#'   constraints_vector = constraints_vector)
#'
#' # Merge the results
#' merge_results(
#'   enriched_pathways, enriched_pathways_directional,
#'   gmt_file = fname_GMT2,
#'   output_prefix = "Aggregated",
#'   col_colors = c("#FF0000", "#00FF00", "#FFFFF0"),
#'   tests = c('rna', 'protein', 'combined')
#' )
#' }
merge_results <- function(enriched_pathways, enriched_pathways_directional, gmt_file, output_prefix = "Aggregated",
                          col_colors = NULL,
                          tests = c(gsub("^Genes_", "", grep("^Genes_", colnames(enriched_pathways), value = TRUE)), 'combined') 
                         ) {
  # read the gmt file
  gmt <- read.GMT(gmt_file)

  if (is.null(tests) || length(tests) == 0) {
    stop("Tests parameter must be provided (e.g., c('rna', 'protein', 'combined')) or NULL")
  }
  
  # check if gmt is a valid GMT object
  if (!is.GMT(gmt)) stop('gmt is not a valid GMT object')

  # Check if all tests (except 'combined') exist in both dataframes with "Genes_" prefix
  tests_to_check <- tests[tests != "combined"]
  if (!all(sapply(tests_to_check, function(test) paste0("Genes_", test) %in% colnames(enriched_pathways))) ||
      !all(sapply(tests_to_check, function(test) paste0("Genes_", test) %in% colnames(enriched_pathways_directional)))) {
    stop("All tests (except 'combined') must exist as columns with 'Genes_' prefix in both enriched_pathways and enriched_pathways_directional")
  }
  
  if (is.null(col_colors)) {
    col_colors <- grDevices::rainbow(length(tests))
  } else if (length(col_colors) != length(tests)) {
    stop("col_colors must have the same length as tests")
  }
  
# Write the combined pathways txt file
df_txtpathways <- data.frame(term_id = c(enriched_pathways$term_id, enriched_pathways_directional$term_id),
                              term_name = c(enriched_pathways$term_name, enriched_pathways_directional$term_name),
                              adjusted_p_val = c(enriched_pathways$adjusted_p_val, enriched_pathways_directional$adjusted_p_val))

# Add the evidence column
df_txtpathways$evidence <- append(enriched_pathways$evidence, enriched_pathways_directional$evidence)

# Aggregate the adjusted p-values by term ID and name
pathways_txt <- stats::aggregate(x = adjusted_p_val ~ term_id + term_name, 
                                  data = df_txtpathways, 
                                  FUN = function(x) {
                                    c(min = min(x))
                                  })

# Convert to data.frame since write.table works natively with data.frames
pathways_txt <- as.data.frame(pathways_txt)
utils::write.table(pathways_txt, 
                   file = paste0(output_prefix, "_pathways.txt"), 
                   row.names = FALSE, 
                   sep = "\t", 
                   quote = FALSE)

# Get unique term IDs and evidence from both analyses
term_ids <- unique(df_txtpathways$term_id)
evidence_combined <- df_txtpathways$evidence

# Create evidence matrix directly
evidence_matrix <- sapply(term_ids, function(id) {
  # Get the evidence list for this pathway
  evidence_list <- evidence_combined[which(df_txtpathways$term_id == id)]
  # Split evidence string if it contains multiple values
  evidence_values <- unlist(strsplit(evidence_list[[1]], ", "))
  # Convert evidence list to binary indicators (0/1) showing which tests are present
  as.numeric(tests %in% evidence_values)
})

# Create final data frame
col_significance <- data.frame(
  term_id = term_ids, 
  t(evidence_matrix), 
  stringsAsFactors = FALSE,
  row.names = NULL
)
names(col_significance) <- c("term_id", tests)

# Check for lost, gained, and shared pathways between methods
lostp <- enriched_pathways$term_id[!enriched_pathways$term_id %in% enriched_pathways_directional$term_id]
gainedp <- enriched_pathways_directional$term_id[!enriched_pathways_directional$term_id %in% enriched_pathways$term_id]
sharedp <- enriched_pathways$term_id[enriched_pathways$term_id %in% enriched_pathways_directional$term_id]

col_significance$directional_impact <- 0
col_significance[col_significance$term_id %in% lostp,]$directional_impact <- 1
col_significance[col_significance$term_id %in% gainedp,]$directional_impact <- 2
  
col_significance <- as.data.table(col_significance)
  
# Add the instruct string for cytoscape
instruct_str <- paste('piechart:',
                      ' attributelist="', 
                      paste(tests, collapse = ','),
                      '" colorlist="', 
                      paste(col_colors, collapse = ','), 
                      '" showlabels=FALSE', sep = '')
col_significance[, "instruct" := instruct_str]
  
utils::write.table(col_significance, 
                   file = paste0(output_prefix, "_subgroups.txt"), 
                   row.names = FALSE, 
                   sep = "\t", 
                   quote = FALSE)

# Filter to include only the specified term IDs
filtered_gmt <- gmt[term_ids]

write.GMT(filtered_gmt, paste0(output_prefix, "_pathways.gmt"))

  return(pathways_txt)
}

#' Merge and filter GMT files based on a list of term IDs
#'
#' This function reads a GMT file and filters it to include only the specified term IDs.
#' It's useful for creating a filtered GMT file for visualization in Cytoscape.
#'
#' @param gmt_file Path to the GMT file
#' @param term_ids Character vector of term IDs to include
#'
#' @return A GMT object containing only the specified terms
#' @export
#'
#' @examples
#' \dontrun{
#' # Get term IDs from merged results
#' merged_results <- merge_results(
#'   enriched_pathways, enriched_pathways_directional,
#'   gmt_file = fname_GMT2,
#'   output_prefix = "Aggregated",
#'   tests = c('rna', 'protein', 'combined'),
#'   col_colors = c("#FF0000", "#00FF00", "#FFFFF0")
#' )
#'
#' # Merge and filter GMT file
#' merge_gmt(
#'   gmt_file = fname_GMT2,
#'   term_ids = merged_results$term_ids
#' )
#' }
merge_gmt <- function(gmt_file, term_ids) {
  # Read the GMT file
  gmt_main <- read.GMT(gmt_file)
  
  # Filter to include only the specified term IDs
  gmt_filtered <- gmt_main[term_ids]
  
  return(gmt_filtered)
}
