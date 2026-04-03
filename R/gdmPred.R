#' calculation of stability using a generalized dissimilarity model.
#'
#' This function will take the community dissimilarity matrix and the
#' environmental metadata table/matrix as input, and make predictions based on
#' a generalized dissimilarity model (GDM), with optional geographic
#' information (X and Y can be longitude and latitude). This model considers
#' the nonlinear relationship between community dissimilarity and environmental
#' distance, and can also include geographic distance as a predictor.
#'
#' @param comdist The community dissimilarity matrix
#' @param envmeta The environmental metadata table/matrix
#' @param sitenames The names of the site (optional, default: NULL)
#' @param X The X coordinates of the sites (optional, default: NULL)
#' @param Y The Y coordinates of the sites (optional, default: NULL)
#' @param geo_enabled Whether to include geographic info  (default: TRUE)
#'
#' @importFrom usedist dist_get
#' @importFrom gdm formatsitepair gdm
#' @returns a column vector of predicted stability values for each site
#'
#' @examples
#' library(vegan)
#' data(varespec)
#' data(varechem)
#' data(BCI)
#' data(BCI.env)
#' example.comdist <- vegdist(varespec)
#' example.stability_GDM <- gdmPred(example.comdist, varechem)
#'
#' example.stability_GDM_geo <- gdmPred(vegdist(BCI, "bray"),
#'     BCI.env[, c("Precipitation", "Elevation", "EnvHet")],
#'     X = BCI.env$UTM.EW,
#'     Y = BCI.env$UTM.NS
#' )
#'
#' @export
gdmPred <- function(
    comdist,
    envmeta,
    sitenames = NULL,
    X = NULL,
    Y = NULL,
    geo_enabled = TRUE
) {
    if (is.null(sitenames)) {
        if (identical(labels(comdist), rownames(envmeta))) {
            sitenames <- labels(comdist)
        } else {
            stop("The labels(comdist) and rownames(envmeta) are not identical!")
        }
    }
    siteids <- if (identical(
        as.character(as.numeric(sitenames)),
        sitenames
    )) {
        as.numeric(sitenames)
    } else {
        seq_len(length(sitenames))
    }

    if (is.null(X) || is.null(Y)) {
        X <- rep(0, length(sitenames))
        Y <- rep(0, length(sitenames))
        geo_enabled <- FALSE
    }

    env_data <- cbind(
        site = siteids,
        X = X,
        Y = Y,
        envmeta
    )
    gdmDissim <- data.frame(
        site = siteids,
        as.matrix(comdist),
        stringsAsFactors = FALSE
    )
    gdm_data <- gdm::formatsitepair(
        bioData = gdmDissim,
        bioFormat = 3,
        predData = env_data,
        siteColumn = "site",
        XColumn = "X",
        YColumn = "Y"
    )
    gdm_model <- gdm(gdm_data, geo = geo_enabled)
    gdm_pred <- predict(gdm_model, data = gdm_data)

    gdm_matrix <- pred2matrix(
        pred = gdm_pred,
        site_ids = sitenames
    )

    result <- data.frame(matrix(NA, nrow = length(labels(comdist)), ncol = 1))
    for (n.site in seq_len(length(sitenames))) {
        sitename <- sitenames[n.site]
        predicted.dist <- mean(gdm_matrix[sitename, ])
        othersites <- setdiff(sitenames, sitename)
        selected.dist <- dist_get(comdist, sitename, othersites)
        mean.measured.dist <- mean(selected.dist)
        result[n.site, 1] <- calcStability(predicted.dist, mean.measured.dist)
    }

    colnames(result)[1] <- "stability_GDM"
    rownames(result) <- sitenames
    return(result)
}
