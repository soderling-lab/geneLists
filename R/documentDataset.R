#' documentDataset
#'
#' generate GMT dataset documentation and save gmt object as .rda
#'
#' @param none
#'
#' @return
#'
#' @author Tyler W Bradshaw, \email{twesleyb10@gmail.com}
#'
#' @references none
#'
#' @keywords
#'
#' @export
#'
#' @examples
#' documentDataset(gmt_file, short_name, Rdir, datadir)
documentDataset <- function(gmt_file, short_name, Rdir, datadir,
                            tags = NULL, save_version = 2) {
  # Add file extension if necessary.
  if (!grepl("\\.gmt", gmt_file)) {
    gmt_file <- paste0(gmt_file, ".gmt")
  }
  gmt_list <- read_gmt(gmt_file, output = "gmt_list")
  # Check if gmt data contains multiple sources.
  n <- length(unique(gmt_list[[2]]))
  if (n > 1) {
    stop("More than one gmt source detected.")
  }
  # Generate documentation.
  dataset <- readLines(file.path(Rdir, "data_emptyDataset.R"))
  namen <- paste(unique(gmt_list$names), collapse = "|")
  dataset[1] <- gsub("dataset_name", namen, dataset[1])
  dataset[13] <- gsub("url", unique(gmt_list$source), dataset[13])
  dataset[17] <- gsub('# \"short_name\"', paste0('"', short_name, '"'), dataset[17])
  myfile <- file.path(Rdir, paste0(short_name, ".R"))
  writeLines(dataset, myfile)
  # Save rda object.
  assign(short_name, gmt_list$genes)
  save(
    list = short_name,
    file = file.path(datadir, paste0(short_name, ".rda")),
    version = save_version
  )
}
