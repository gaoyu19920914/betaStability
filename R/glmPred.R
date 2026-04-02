#' calculation of stability using a generalized linear model.
#'
#' This function will take the dissimilarity matrix and the environmental
#' matrix as input, and calculate the stability of each site using a generalized
#' linear model (gLM), where the contributions are constrained as non-negative
#'  `lower.limits=0` to ensure the explanability of each coefficient.
#' The stability is calculated by comparing the predicted distance (based on the
#' linear model) and the mean measured distance (based on vegdist function).
#'
#' @param comdist The community dissimilarity matrix
#' @param envmeta The environmental metadata table/matrix
#' @param sitenames The names of the site
#'
#' @importFrom usedist dist_subset dist_get
#' @importFrom BBmisc normalize
#' @importFrom glmnet cv.glmnet
#' @returns a column vector of predicted stability values for each site
#'
#' @examples
#' data(varespec)
#' data(varechem)
#' example.comdist <- vegdist(varespec)
#' example.stability_GLM <- glmPred(example.comdist, varechem)
#'
#' @export
glmPred <- function(comdist,
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
      # print(i)
      # print(j)
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
    glm_predictors <- subset(trainingset, select = names(deltaenvnorm))
    glm_beta <- subset(trainingset, select = c("beta"))
    this.cv.glmnet <- cv.glmnet(as.matrix(glm_predictors),
                                as.matrix(glm_beta),
                                lower.limits = 0)
    best_lambda <- this.cv.glmnet$lambda.min
    coef(this.cv.glmnet, s = best_lambda)
    beta_pred <- predict(this.cv.glmnet,
                         newx = as.matrix(subset(validatingset,
                                                 select=names(deltaenvnorm))),
                         s = best_lambda)
    # plot(as.matrix(subset(validatingset, select = c("beta"))),
    #      beta_pred,
    #      col = "blue",
    #      xlab = "Observed Beta Diversity Indices",
    #      ylab = "GLM Predicted Beta Diversity Indices")
    # abline(0, 1, col = "red")
    othersites <- setdiff(sitenames, sitename)
    selected.dist <- dist_get(comdist, sitename, othersites)
    mean.measured.dist <- mean(selected.dist)
    result[n.site, 1] <- calcStability(mean(beta_pred), mean.measured.dist)
  }
  colnames(result)[1] <- "stability_GLM"
  rownames(result) <- sitenames
  return(result)
}
