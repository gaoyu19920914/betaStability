#' calculation of stability using a xgboost model.
#'
#' This function will take the dissimilarity matrix and the environmental
#' matrix as input, and calculate the stability of each site using a xgboost
#' model to improve the prediction performance.
#'
#' @param comdist The community dissimilarity matrix
#' @param envmeta The environmental metadata table/matrix
#' @param sitenames The names of the site
#' @param seed The random seed for reproducibility of the xgboost model
#' @param params A list of parameters for the xgboost model.
#' If NULL, default parameters will be used.
#'
#' @importFrom usedist dist_subset dist_get
#' @importFrom BBmisc normalize
#' @importFrom xgboost xgb.DMatrix xgb.train
#' @returns a column vector of predicted stability values for each site
#'
#' @examples
#' data(varespec)
#' data(varechem)
#' example.comdist <- vegdist(varespec)
#' example.stability_XGB <- xgboostPred(example.comdist, varechem)
#'
#' @export
xgboostPred <- function(comdist,
    envmeta,
    sitenames = NULL,
    seed = NULL,
    params = NULL) {
    if (!is.null(seed)) set.seed(seed)
    if (is.null(params)) {
        params <- list(
            objective = "reg:squarederror",
            max_depth = 6,
            eta = 0.1,
            subsample = 0.8,
            colsample_bytree = 0.8,
            gamma = 0,
            min_child_weight = 1
        )
    }

    result <- data.frame(matrix(NA, nrow = length(labels(comdist)), ncol = 1))
    if (is.null(sitenames)) {
        if (identical(labels(comdist), rownames(envmeta))) {
            sitenames <- labels(comdist)
        } else {
            stop("The labels(comdist) and rownames(envmeta) are not identical!")
        }
    }

    # prepare df.ij.deltaenv.beta dataframe between pairs of sites.
    df.ij.deltaenv.beta <- data.frame(
        matrix(ncol = 2 + ncol(envmeta) + 1, nrow = 0))
    colnames(df.ij.deltaenv.beta) <- c("i", "j", colnames(envmeta), "beta")
    for (i in seq_len(nrow(envmeta) - 1)) {
        for (j in (i + 1):nrow(envmeta)) {
            row.to.add <- unlist(c(
                i,
                j,
                get_deltaenv_rows(i, j, envmeta),
                dist_get(comdist, i, j)
            ))
            row.to.add <- setNames(row.to.add, colnames(df.ij.deltaenv.beta))
            df.ij.deltaenv.beta <- rbind(
                df.ij.deltaenv.beta,
                row.to.add
            )
        }
    }
    colnames(df.ij.deltaenv.beta) <- c("i", "j", colnames(envmeta), "beta")

    deltaenvnorm <- normalize(
        df.ij.deltaenv.beta[, !names(df.ij.deltaenv.beta) %in%
            c("i", "j", "beta")],
        method = "range",
        margin = 2
    )

    df.ij.deltaenv.beta.norm <- cbind(
        df.ij.deltaenv.beta[, c("i", "j")],
        deltaenvnorm,
        subset(df.ij.deltaenv.beta,
            select = c("beta")
        )
    )

    for (n.site in seq_len(nrow(envmeta))) {
        sitename <- sitenames[n.site]
        validatingset <- df.ij.deltaenv.beta.norm[
            df.ij.deltaenv.beta.norm$i == n.site |
                df.ij.deltaenv.beta.norm$j == n.site,
        ]
        trainingset <- df.ij.deltaenv.beta.norm[
            !(df.ij.deltaenv.beta.norm$i == n.site |
                df.ij.deltaenv.beta.norm$j == n.site),
        ]
        xgb_beta <- trainingset$beta
        xgb_features <- trainingset[, !names(trainingset) %in%
            c("i", "j", "beta")]
        validatingset_features <- validatingset[, !names(validatingset) %in%
            c("i", "j", "beta")]

        # Convert data to xgb.DMatrix format
        dtrain <- xgboost::xgb.DMatrix(
            data = as.matrix(xgb_features),
            label = xgb_beta
        )
        dtest <- xgboost::xgb.DMatrix(data = as.matrix(validatingset_features))

        # Train xgboost model
        model_xgb <- xgb.train(
            params = params,
            data = dtrain,
            nrounds = 100,
            verbose = 0,
            evals = list(train = dtrain),
            print_every_n = 50
        )

        beta_pred <- predict(model_xgb, dtest)

        othersites <- setdiff(sitenames, sitename)
        selected.dist <- dist_get(comdist, sitename, othersites)
        mean.measured.dist <- mean(selected.dist)
        result[n.site, 1] <- calcStability(mean(beta_pred), mean.measured.dist)
    }
    colnames(result)[1] <- "stability_XGB"
    rownames(result) <- sitenames
    return(result)
}
