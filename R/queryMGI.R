#' queryMGI
#'
#' Query a saved copy of the MGI gene database for gene identifiers. 
#' MGI gene_map includes MGI, UniProt Accession, and Entrez identifiers.
#' Maps UniProt Accession to Entrez by default. Arguments 'from' and 'to' are
#' case-insensitive.
#'
#' See also `data(mgi_map)` # MGI | UniProt | Entrez mapping df
#'
#' @author Tyler W Bradshaw, \email{twesleyb10@gmail.com}
#'
#' @references http://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt
#'
#' @keywords MGI batch query
#'
#' @export queryMGI

queryMGI <- function(input_ids, from = "UniProt", to = "Entrez") {

  # map gene identifiers using MGI mapping data

  # load stored copy of mgi gene mapping data
  data(list="mgi_map", package="geneLists") # mgi_map

  gene_map <- mgi_map

  # get column cooresponding to input gene identifier type (from)
  idy <- grep(tolower(from), tolower(colnames(gene_map)))
  if (!is.numeric(idy) | length(idy)!=1) {
	  stop("Input identifiers must be one of: ", paste(colnames(gene_map),sep=", "))
  }
  id_from <- colnames(gene_map)[idy]

  idy <- grep(tolower(to), tolower(colnames(gene_map)))
  if (!is.numeric(idy) | length(idy)!=1) {
	  stop("Output identifiers must be one of: ", paste(colnames(gene_map),sep=", "))
  }
  id_to <- colnames(gene_map)[idy]


  # map IDs
  idx <- match(input_ids, gene_map[[id_from]])

  output_ids <- gene_map[[id_to]][idx]

  # check for missing
  n_missing <- sum(is.na(output_ids))
  if (n_missing > 0) {
    warning("Unable to map ", n_missing, " gene identifiers to entrez.")
  }

  # return input identifiers mapped to entrez
  names(output_ids) <- input_ids

  return(output_ids)
} #EOF
