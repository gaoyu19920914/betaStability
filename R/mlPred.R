#' calculation of stability using multiple linear regression model.
#'
#' This function will take the diversity matrix and the environmental distance
#' matrix as input, and calculate the stability of each site using multiple
#' linear model (ML). The stability is calculated by comparing the predicted
#' distance (based on the multiple linear model) and the mean measured distance
#' between the site and other sites (based on the difference of envmeta and
#' the corresponding comdist).
#'
#' @param comdist The community dissimilarity matrix
#' @param envmeta The environmental metadata table/matrix
#' @param sitenames The names of the site
#'
#' @importFrom usedist dist_subset dist_get
#' @importFrom BBmisc normalize
#' @returns a column vector of predicted stability values for each site
#'
#' @examples
#' data(varespec)
#' data(varechem)
#' example.comdist <- betaDiv(varespec)
#' example.stability_ML <- mlPred(example.comdist, varechem)
#'
#' @export
mlPred <- function(comdist,
                   envmeta,
                   sitenames=NULL) {
  result <- data.frame(matrix(NA, nrow = length(labels(comdist)), ncol =1))
  if (is.null(sitenames)){
    if(identical(labels(comdist), rownames(envmeta))){
      sitenames <- labels(comdist)
    } else {
      stop("The labels of comdist and rownames of envmeta are not identical!")
    }
  }

  # prepare df.ij.deltaenv.beta dataframe between pairs of sites.
  df.ij.deltaenv.beta <- data.frame(matrix(ncol = 2+ncol(envmeta)+1, nrow = 0))
  colnames(df.ij.deltaenv.beta) <- c("i", "j", colnames(envmeta), "beta")

  for (i in 1:(nrow(envmeta)-1)) {
    for (j in (i+1):nrow(envmeta)){
      row.to.add <- unlist(c(i,
                             j,
                             get_deltaenv_rows(i, j, envmeta),
                             dist_get(comdist, i, j)))
      row.to.add <- setNames(row.to.add, colnames(df.ij.deltaenv.beta))
      df.ij.deltaenv.beta <- rbind(df.ij.deltaenv.beta,
                                   row.to.add)
    }
  }

  colnames(df.ij.deltaenv.beta) <- c("i", "j", colnames(envmeta), "beta")
  deltaenvnorm <- normalize(
    df.ij.deltaenv.beta[, !names(df.ij.deltaenv.beta) %in% c("i", "j", "beta")],
    method = "range",
    margin = 2)

  df.ij.deltaenv.beta.norm <- cbind(df.ij.deltaenv.beta[, c("i", "j")],
                                    deltaenvnorm,
                                    subset(df.ij.deltaenv.beta,
                                           select = c("beta")))
  for (n.site in 1:(nrow(envmeta))) {
    sitename <- sitenames[n.site]
    validatingset <- df.ij.deltaenv.beta.norm[
      df.ij.deltaenv.beta.norm$i == n.site |
        df.ij.deltaenv.beta.norm$j == n.site ,]
    trainingset <- df.ij.deltaenv.beta.norm[
      !(df.ij.deltaenv.beta.norm$i == n.site |
          df.ij.deltaenv.beta.norm$j == n.site) ,]
    mlm_variables <- paste(names(deltaenvnorm), collapse = " + ")
    mlm_formula <- as.formula(paste("beta ~", mlm_variables))
    this.mlm <- lm(mlm_formula, data = trainingset)
    predicted.dist <- mean(predict(this.mlm, validatingset))
    othersites <- setdiff(sitenames, sitename)
    selected.dist <- dist_get(comdist, sitename, othersites)
    mean.measured.dist <- mean(selected.dist)
    result[n.site, 1] <- calcStability(predicted.dist, mean.measured.dist)
  }
  colnames(result)[1] <- "stability_ML"
  rownames(result) <- sitenames
  return(result)
}
