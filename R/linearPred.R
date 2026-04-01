#' calculation of stability using linear prediction model.
#'
#' This function will take the diversity matrix and the environmental distance
#' matrix as input, and calculate the stability of each site using linear model.
#' The stability is calculated by comparing the predicted distance (based on the
#' linear model) and the mean measured distance (based on betaDiv function).
#'
#' @param comdist The community dissimilarity matrix
#' @param envdist The environmental dissimilarity matrix
#' @param sitenames The names of the site
#'
#' @importFrom usedist dist_subset dist_get
#' @returns a column vector of predicted stability values for each site
#'
#' @examples
#' data(varespec)
#' data(varechem)
#' example.comdist <- betaDiv(varespec)
#' example.envdist <- dist(
#'   BBmisc::normalize(
#'     varechem,
#'     method = "range",
#'     margin = 2),
#'   method = "euclidean")
#' example.stability_LM <- linearPred(example.comdist, example.envdist)
#'
#' @export
linearPred <- function(comdist,
                       envdist,
                       sitenames=NULL) {

  result <- data.frame(matrix(NA, nrow = length(labels(comdist)), ncol =1))

  if (is.null(sitenames)){
    if(identical(labels(comdist), labels(envdist))){
      sitenames <- labels(comdist)
    } else {
      stop("The labels of comdist and envdist are not identical!")
    }
  }

  for (n.site in 1:length(sitenames)){
    sitename <- sitenames[n.site]
    subcomdist <- dist_subset(comdist, labels(comdist) != sitename)
    subenvdist <- dist_subset(envdist, labels(envdist) != sitename)

    y = unlist(as.list(subcomdist))
    x = unlist(as.list(subenvdist))

    this.linear.model <- lm(y ~ x)

    selected.comdist <- dist_get(comdist, sitename, sitenames)
    mean.dist <- mean(selected.comdist)
    selected.envdist <- dist_get(envdist, sitename, sitenames)
    mean.envdist <- mean(selected.envdist)

    predicted.dist <- predict(this.linear.model,
                              newdata = data.frame(x = mean.envdist))
    result[n.site,1] <- calcStability(predicted.dist, mean.dist)
  }

  colnames(result)[1] <- "stability_Linear"
  rownames(result) <- sitenames
  return(result)
}
