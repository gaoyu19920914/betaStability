#   ##### 3.4 gdm nonlinear prediction ####
#' calculation of stability using a (nonlinear) generalized dissimilarity model.
#'
#' TODO: fill this description
#' This function ......
#'
#' @param comdist The community dissimilarity matrix
#' @param envmeta The environmental metadata table/matrix
#' @param sitenames The names of the site (optional, default: NULL)
#' @param X The X coordinates of the sites (optional, default: NULL)
#' @param Y The Y coordinates of the sites (optional, default: NULL)
#'
#' @importFrom usedist dist_get
#' @importFrom gdm formatsitepair gdm
#' @returns a column vector of predicted stability values for each site
#'
#' @examples
#' varespec <- data(varespec)
#' varechem <- data(varechem)
#' example.comdist <- betaDiv(varespec)
#' example.stability_GDM <- gdmPred(example.comdist, varechem)
#'
#' @export
gdmPred <- function(comdist,
                                  envmeta,
                                  sitenames=NULL,
                                  X=NULL,
                                  Y=NULL) {
  if (is.null(sitenames)){
    if(identical(labels(comdist), rownames(envmeta))){
      sitenames <- labels(comdist)
    } else {
      stop("The labels of comdist and rownames of envmeta are not identical!")
    }
  }
  siteids <- if (identical(as.character(as.numeric(sitenames)),
                           sitenames)){
    as.numeric(sitenames)
  } else {
    1:length(sitenames)
  }

  if(is.null(X) || is.null(Y)){
    X <- rep(0, length(sitenames))
    Y <- rep(0, length(sitenames))
  }

  env_data <- cbind(site = siteids,
                    X = X,
                    Y = Y,
                    envmeta)
  gdmDissim <- data.frame(site = siteids,
                          as.matrix(comdist),
                          stringsAsFactors = FALSE)
  gdm_data <- formatsitepair(
    bioData = gdmDissim,
    bioFormat = 3,
    predData = env_data,
    siteColumn = "site",
    XColumn = "X",
    YColumn = "Y"
  )
  gdm_model <- gdm(gdm_data, geo = FALSE)
  gdm_pred <- predict(gdm_model, data = gdm_data)

  gdm_matrix <- pred2matrix(pred = gdm_pred,
                            site_ids = sitenames)

  result <- data.frame(matrix(NA, nrow = length(labels(comdist)), ncol =1))
  for (n.site in 1:length(sitenames)) {
    sitename <- sitenames[n.site]
    predicted.dist <- mean(gdm_matrix[sitename,])
    othersites <- setdiff(sitenames, sitename)
    selected.dist <- dist_get(comdist, sitename, othersites)
    mean.measured.dist <- mean(selected.dist)
    result[n.site, 1] <- calcStability(predicted.dist, mean.measured.dist)
  }

  colnames(result)[1] <- "stability_GDM"
  rownames(result) <- sitenames
  return(result)
}
