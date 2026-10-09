#' calculate stability based on predicted and measured distances.
#'
#' This function calculates the stability of a site by comparing the predicted
#' distance and the measured distance.
#' The algorithm is simple and straightforward: if the measured distance is
#' greater than the predicted distance, the stability is negative, and if the
#' measured distance is less than the predicted distance, the stability is
#' positive. The stability value is normalized to be between -1 and 1, where -1
#' indicates the least stable (measured distance is much greater than predicted)
#' and 1 indicates the most stable (measured distance is much less than
#' predicted).
#' To avoid magnifying the deviations near 0 and 1 of the predicted.dist, a
#' symmetric algorithm is introduced without changing the range of the
#' calculated stability value (-1 to 1).
#'
#' @param predicted.dist The predicted distance (shall be in range 0~1)
#' @param measured.dist The measured distance (shall be in range 0~1)
#' @param symmetric Whether to use symmetric algorithm in the calculating
#' (optional, default: FALSE).
#'
#' @returns a numeric value of stability for the site in range \[-1, 1\].
#'
#' @examples
#' calcStability(predicted.dist = 0.3, measured.dist = 0.5)
#' calcStability(predicted.dist = 0.3, measured.dist = 0.1)
#' calcStability(predicted.dist = 0.3, measured.dist = 0.5, symmetric = TRUE)
#' calcStability(predicted.dist = 0.3, measured.dist = 0.1, symmetric = TRUE)
#'
#' @export
calcStability <- function(predicted.dist, measured.dist, symmetric=FALSE) {
  if (!is.numeric(predicted.dist) || !is.numeric(measured.dist)) {
    stop("Error: predicted.dist and measured.dist should be numeric.")
  }
  if (predicted.dist < 0 || predicted.dist > 1) {
    stop("Error: predicted.dist should be in range [0, 1]")
  }
  if (measured.dist < 0 || measured.dist > 1) {
    stop("Error: measured.dist should be in range [0, 1]")
  }
  if(symmetric){
    return ((predicted.dist - measured.dist) / (predicted.dist + measured.dist))
  } else {
    if (measured.dist > predicted.dist) {
      return(-(measured.dist - predicted.dist) / (1 - predicted.dist))
    } else {
      return((predicted.dist - measured.dist) / predicted.dist)
    }
  }
}
