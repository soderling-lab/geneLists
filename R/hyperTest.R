#' hyperTest
#'
#' Perform hypergeometric test for enrichment
#'
#' @param listA - vector of items to test for enrichment in listB
#'
#' @param listB - vector of items
#'
#' @return p value and fold enrichment
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
#' hyperTest(listA, listB, background)
hyperTest <- function(listA, listB, background, enrichment = TRUE) {
  # hyperTest - perform hypergeometric test.
  listA <- listA[listA %in% background]
  listB <- listB[listB %in% background]
  n_success <- length(intersect(listA, listB))
  N_success <- length(intersect(listA, background))
  N <- length(background)
  n <- length(listB)
  FE <- n_success / (n * (N_success / N))
  if (enrichment) {
    # Test for enrichment
    p <- phyper(n_success - 1, N_success, N - N_success,
      n,
      lower.tail = FALSE
    )
    return(c("Fold enrichment" = FE, "P-value" = p))
  } else {
    # Test for depletion.
    p <- phyper(n_success, N_success, N - N_success,
      n,
      lower.tail = TRUE
    )
    return(c("Fold depletion" = 1 / FE, "P-value" = p))
  }
}
