#' mapIDs
#'
#' A simple function that maps gene identifiers using a given gene_map 
#' mapping data.frame. Defaults to uniprot_map. Column names must be typed to
#' match exactly.
#'
#' @export mapIDs

mapIDs <- function(ids, from, to, gene_map=NULL) {

	if (is.null(gene_map)) {
	  data(list="uniprot_map",package="geneLists")
	  gene_map <- uniprot_map
	}

	idx <- match(ids,gene_map[[from]])

	return(gene_map[[to]][idx])
} #EOF
