#' read_gmt - read gmt file into R
#'
#' returns a list of pathways from a GMT file.
#'
#' @param gmt_file Path to a GMT file.
#'
#' @return A list of vectors with gene sets.
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords
#'
#' @importFrom  utils head tail
#'
#' @export
#'
#' @examples
#' read_gmt("path2/gmt_file")
read_gmt <- function(gmt_file, output = "gene_list") {
  gmt_lines <- strsplit(readLines(gmt_file), "\t")
  gene_list <- lapply(gmt_lines, tail, -2)
  gmt_names <- sapply(gmt_lines, head, 1)
  names(gene_list) <- gmt_names
  gmt_source <- sapply(gmt_lines, "[", 2)
  gmt_list <- list("names" = gmt_names, "source" = gmt_source, "genes" = gene_list)
  if (output == "gene_list") {
    return(gene_list)
  } else if (output == "gmt_list") {
    return(gmt_list)
  }
}
