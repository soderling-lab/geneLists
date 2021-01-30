#' geneLists
#'
#' list all available datasets, search for a dataset given a tag.
#'
#' @param none
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
#' geneLists()
geneLists <- function(pattern = NULL) {
  # List all available geneLists.
  gene_lists <- data(package = "geneLists")$results[, "Item"]
  names(gene_lists) <- NULL
  if (!is.null(pattern)) {
    # If provided, search for user specified pattern.
    return(gene_lists[grep(pattern, gene_lists)])
  } else {
    # Otherwise, return all available datastes.
    return(gene_lists)
  }
}
