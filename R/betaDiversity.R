#' calculation of beta diversity between all pairs of communities.
#'
#' This is the first step for all calculations. Beta diversity will be
#' calculated using the vegan::vegdist function.
#'
#' @param x The community matrix
#' @param method The method to calculate beta diversity
#' @param binary Whether convert to presence/absence (using decostand) before
#' calculation
#' @param diag Compute diagonals
#' @param upper Whether return only the upper diagonal in the matrix
#' @param ... Other parameters to parse to vegan::vegdist()
#'
#' @importFrom vegan vegdist
#' @returns a dissimilarity matrix
#'
#' @examples
#' ...
#'
#' @export
betaDiv <- function(x,
                    method = "bray",
                    binary = FALSE,
                    diag = FALSE,
                    upper = FALSE,
                    ...) {
  result <- vegdist(x, method, binary, diag, upper, ...)
  return(result)
}

