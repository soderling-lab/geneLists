#' write_gmt
#'
#' write list to gmt file
#'
#' @param gmt_list - named list of entrez ids
#'
#' @return
#'
#' @author Tyler W Bradshaw, \email{tyler.w.bradshaw@duke.edu}
#'
#' @references none
#'
#' @keywords
#'
#' @export
#'
#' @examples
#' write_gmt(gmt_list, gmt_source, gmt_file)
write_gmt <- function(gmt_list, gmt_source, gmt_file) {
  if (!is.list(gmt_list)) {
    stop("Please provide a named list of genes as input.")
  }
  gmt_names <- names(gmt_list)
  gmt_genes <- sapply(gmt_list, function(x) paste(x, collapse = "\t"))
  gmt <- paste(gmt_names, gmt_source, gmt_genes, sep = "\t")
  # Add gmt extension to file if it isnt present.
  if (!tools::file_ext(gmt_file) == "gmt") {
    gmt_file <- paste0(gmt_file, ".gmt")
  }
  myfile <- file(gmt_file)
  writeLines(gmt, myfile)
  close(myfile)
}
