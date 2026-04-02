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
#'
#' @param predicted.dist The predicted distance (shall be in range 0~1)
#' @param measured.dist The measured distance (shall be in range 0~1)
#'
#' @returns a numeric value of stability for the site in range \[-1, 1\].
#'
#' @examples
#' calcStability(predicted.dist = 0.3, measured.dist = 0.5)
#' calcStability(predicted.dist = 0.3, measured.dist = 0.1)
#'
#' @export
calcStability <- function(predicted.dist, measured.dist) {
    if (!is.numeric(predicted.dist) || !is.numeric(measured.dist)) {
        message("Error: predicted.dist and measured.dist should be numeric.")
        return(NA)
    }
    if (predicted.dist < 0 || predicted.dist > 1) {
        message("Error: predicted.dist should be in range [0, 1]")
        return(NA)
    }
    if (measured.dist < 0 || measured.dist > 1) {
        message("Error: measured.dist should be in range [0, 1]")
        return(NA)
    }
    if (measured.dist > predicted.dist) {
        return(-(measured.dist - predicted.dist) / (1 - predicted.dist))
    } else {
        return((predicted.dist - measured.dist) / predicted.dist)
    }
}
