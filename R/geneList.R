#' geneList
#'
#' a constructor function for the "geneList" class
#'
#' @export geneList

geneList <- function(name, genes, species=NULL, reference=NULL) {
		# checks
		stopifnot("bibentry" %in% class(reference))
		stopifnot(is.vector(genes))
		# assemble gene list
		obj <- list()
		obj[[name]] <- genes
		attr(obj, "class") <- "geneList"
		attr(obj, "species") <- species
		attr(obj, "bibtex") <- reference
		return(obj)
}
