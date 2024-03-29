#' getHomologs
#'
#' Download the NCBI homology database and map input Entrez ids to homologs in
#' a given species of interest.
#'
#' @param hitpredict (data.table) - HitPredict data.
#'
#' @param entrez (character) - Entrez IDs -- can be from any species
#'
#' @return none
#'
#' @author Tyler W Bradshaw, \email{twesleyb10@gmail.com}
#'
#' @references none
#'
#' @keywords none
#'
#' @export getHomologs
#'
#' @importFrom data.table fread 
#' @importFrom dplyr %>% filter select rename_all

#' @examples
#' getHomologs(entrez, species = "mouse")


getHomologs <- function(entrez, species = NULL, taxid = NULL, quiet = TRUE) {

  # homolgene data URL
  URL <- "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"

  # Get taxid of provided species.
  dbs <- mappingDBs()

  # Check that user has provided a species or taxid.
  if (is.null(species) & is.null(taxid)) {
    stop(paste("Please provide either a taxid or species alias."))
  }

  # Check if species is valid.
  if (!is.null(species)) {
    all_species <- sapply(dbs, "[[", "alias")
    idx <- grep(species, sapply(dbs, "[[", "alias"))
    if (length(idx) == 0) {
      stop(paste0(
        "Please provide a valid species taxid or alias: ",
        paste(all_species, collapse = ", "), "."
      ))
    }
    taxid <- dbs[[idx]][["taxid"]]
  }

  # Check if taxid is valid.
  if (!is.null(taxid)) {
    is_valid <- taxid %in% sapply(dbs, "[[", 1)
    if (!is_valid) {
      stop("Please provide a valid taxid.")
    }
  }

  # status
  if (!quiet) {
    message("Downloading NCBI HomoloGene data...")
  }

  # Download and load NCBI homology gene data.
  destfile <- file.path(getwd(), "homologene.data")
  download.file(URL, destfile, quiet = quiet)
  homologene <- data.table::fread(destfile, header = FALSE)

  # Remove raw data.
  unlink(destfile)

  # Fix column names.
  # Gene ID is a genes organism specific Entrez ID.
  # HID is the genes homology id.
  homologene <- homologene %>% dplyr::rename_all(list(~ c(
    "HID", "TaxonomyID", "GeneID",
    "GeneSymbol", "ProteinGI", "ProteinAccession"
  )))

  # Use HomoloGene to create homology map.
  homology_map <- as.list(homologene$HID)
  names(homology_map) <- homologene$GeneID

  # Map Entrez to homology ID.
  hid <- homology_map[entrez]
  hid[seq(1:length(hid))[unlist(lapply(hid, is.null))]] <- NA
  # Subset homologene data, keep all genes from your species of interest.
  homologene <- homologene %>% dplyr::filter(TaxonomyID == taxid)
  # Taxonomy info.
  annotationDBs <- mappingDBs()
  osDB <- unlist(annotationDBs[sapply(annotationDBs, "[", 1) == taxid])
  names(osDB) <- sapply(strsplit(names(osDB), "\\."), "[", 2)
  organism <- osDB["alias"]
  # Get entrez associated with these HIDs.
  osEntrez <- homologene$GeneID[match(hid, homologene$HID)]
  # Status report.
  is_missing <- is.na(osEntrez)
  n_total <- length(is_missing)
  if (!quiet) {
    message(round(100 * sum(!is_missing) / n_total, 2),
      "% of supplied entrez IDs were successfully ",
      "mapped to homologous genes in ", organism, ".")
  }
  # Return data with genes mapped to gene specific entrez.
  return(osEntrez)
}
