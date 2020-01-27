#' read_gmt
#' returns a list of pathways from a GMT file.
#' @param gmt.file Path to a GMT file.
#' @return A list of vectors with gene sets.
#' @export
#' @examples
#' pathways <- gmtPathways(system.file(
#'      "extdata", "mouse.reactome.gmt", package="fgsea"))
#' @importFrom  utils head tail
#' @export
gmtPathways <- function(gmt.file) {
	pathwayLines <- strsplit(readLines(gmt.file), "\t")
	pathways <- lapply(pathwayLines, tail, -2)
	names(pathways) <- sapply(pathwayLines, head, 1)
	return(pathways)
}
