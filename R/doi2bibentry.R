#' doi2bibentry
#'
#' get bibtex entry for given doi using rcrossref
#'
#' @export doi2bibentry
#' @importFrom rcrossref cr_cn

doi2bibentry <- function(doi, gsub_key = TRUE) {
  # use rcrossref to get bibtex
  bib <- rcrossref::cr_cn(dois = doi, format = "bibentry")
  # to coerce to class bibentry, 'entry' should be 'bibtype'
  names(bib)[names(bib) == "entry"] <- "bibtype"
  # I prefer AuthorYear formatted keys
  if (gsub_key) {
    # rm underscore from key
    bib$key <- gsub(pattern = "_", replacement = "", bib$key)
  }
  # coerce to class bibentry
  return(do.call(bibentry, bib))
}
