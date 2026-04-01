#   ##### 3.5 GAMs (mgcv) prediction ####
#' calculation of stability using an generalized additive model.
#'
#' TODO: fill this description
#' This function ......
#'
#' @param comtable The community table
#' @param envmeta The environmental metadata table/matrix
#' @param comdist The community dissimilarity matrix (optional, default: NULL)
#' @param sitenames The names of the site (optional, default: NULL)
#' @param GAM.dist.method The method for calculating dist (default: "manhattan")
#'
#' @importFrom usedist dist_get
#' @importFrom mgcv gam
#' @returns a column vector of predicted stability values for each site
#'
#' @examples
#' varespec <- data(varespec)
#' varechem <- data(varechem)
#' example.stability_GAM <- gamPred(varespec, varechem)
#'
#' @export
gamPred <- function(comtable,
                    envmeta,
                    comdist=NULL,
                    sitenames=NULL,
                    GAM.dist.method = "manhattan") {
  if(is.null(comdist)) {
    comdist <- betaDiv(comtable)
  }
  if (is.null(sitenames)){
    if(identical(labels(comdist), rownames(envmeta))){
      sitenames <- labels(comdist)
    } else {
      stop("The labels of comdist and rownames of envmeta are not identical!")
    }
  }

  site_ids <- rownames(comtable)
  y_GAM <- as.matrix(comdist)[lower.tri(comdist)]
  varnames <- colnames(envmeta)
  x_GAM <- lapply(varnames, function(v) {
    m <- as.matrix(dist(varechem[[v]], method = GAM.dist.method))
    m[lower.tri(m)]
  })
  x_GAM <- do.call(cbind, x_GAM)
  colnames(x_GAM) <- varnames
  data_GAM <- cbind(y = y_GAM, x_GAM)
  # TODO: make variables a parameter in function
  formula_GAM_str <- paste("y ~",
                           paste("s(", varnames, ")", sep="", collapse=" + "))
  # TODO: add different methods for GAM predictors
  model_GAM <- gam(as.formula(formula_GAM_str),
                   data = as.data.frame(data_GAM),
                   method = "REML")
  # view the smooth functions of each parameter in varechem
  # plot(model_GAM, pages = 1)

  # predict Y and compare with existing Y:
  pred_y_GAM <- predict(model_GAM,
                        newdata = as.data.frame(data_GAM),
                        type = "response")
  pred_y_GAM_matrix <- pred2matrix(pred = pred_y_GAM,
                                   site_ids = site_ids)
  result <- data.frame(matrix(NA, nrow = length(labels(comdist)), ncol =1))
  for (n.site in 1:length(sitenames)) {
    sitename <- sitenames[n.site]
    predicted.dist <- mean(pred_y_GAM_matrix[sitename,])
    othersites <- setdiff(sitenames, sitename)
    selected.dist <- dist_get(comdist, sitename, othersites)
    mean.measured.dist <- mean(selected.dist)
    result[n.site, 1] <- calcStability(predicted.dist, mean.measured.dist)
  }

  colnames(result)[1] <- "stability_GAM"
  rownames(result) <- sitenames
  return(result)
}
