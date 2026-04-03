#' Beta Stability Calculation with Multiple Prediction Methods
#'
#' This function integrates various prediction methods (linear, multiple linear,
#' generalized linear model, generalized additive model, generalized
#' dissimilarity model, random forest, and xgboost) to calculate the stability
#' of each site.
#'
#' @param comtable The community table (required).
#' @param envmeta The environmental metadata table/matrix (required).
#' @param comdist The community dissimilarity matrix (optional). If not
#' provided, computed from comtable using vegdist().
#' @param envdist The environmental dissimilarity matrix (optional). If not
#' provided, computed from envmeta using euclidean distance on range-normalized
#' envmeta.
#' @param sitenames The names of the site (optional, default: NULL)
#' @param method A character string or vector specifying the prediction
#' method(s) to use. Available options: "linearPred", "mlPred", "glmPred",
#' "gamPred", "gdmPred", "rfPred", "xgboostPred". Use "all" to run all methods.
#' Default is "linearPred".
#' @param X The X coordinates of the sites for gdmPred (optional, default: NULL)
#' @param Y The Y coordinates of the sites for gdmPred (optional, default: NULL)
#' @param geo_enabled Whether to include geographic info for gdmPred
#' (default: TRUE)
#' @param GAM.dist.method The method for calculating dist for gamPred
#' (default: "manhattan")
#' @param vegdist.method The method for calculating comdist from comtable
#' (default: "bray")
#' @param seed The random seed for reproducibility of rfPred and xgboostPred
#' (default: 42)
#' @param xgboost.params A list of parameters for the xgboost model
#' (default: NULL). If NULL, default parameters will be used.
#'
#' @importFrom stats dist
#' @importFrom vegan vegdist
#'
#' @returns If method = "all", returns a data frame with 7 columns, each
#' representing results from each selected method. If method length is 1,
#' returns a column vector of predicted stability values. If method length > 1,
#' returns a data frame with each column representing results from each
#' selected method.
#'
#' @examples
#' library(vegan)
#' data(varespec)
#' data(varechem)
#'
#' # Single method (linearPred)
#' result_linear <- betaStability(
#'     comtable = varespec,
#'     envmeta = varechem,
#'     method = "linearPred"
#' )
#'
#' # Multiple methods
#' results_multi <- betaStability(
#'     comtable = varespec,
#'     envmeta = varechem,
#'     method = c("linearPred", "mlPred")
#' )
#'
#' @export
betaStability <- function(comtable = NULL,
    envmeta = NULL,
    comdist = NULL,
    envdist = NULL,
    sitenames = NULL,
    method = "linearPred",
    X = NULL,
    Y = NULL,
    geo_enabled = TRUE,
    GAM.dist.method = "manhattan",
    vegdist.method = "bray",
    xgboost.params = NULL,
    seed = 42) {
    # Validate required parameters
    if (is.null(comtable)) {
        stop("'comtable' is required")
    }
    if (is.null(envmeta)) {
        stop("'envmeta' is required")
    }

    # Validate method parameter
    valid_methods <- c(
        "linearPred", "mlPred", "glmPred", "gamPred",
        "gdmPred", "rfPred", "xgboostPred"
    )

    # Handle "all" case
    if (length(method) == 1 && method == "all") {
        method <- valid_methods
    } else if (!all(method %in% valid_methods)) {
        invalid_methods <- setdiff(method, valid_methods)
        stop(
            "Invalid method(s): ", paste(invalid_methods, collapse = ", "),
            ". Valid methods are: ", paste(valid_methods, collapse = ", ")
        )
    }

    # Auto-compute comdist from comtable if not provided
    if (is.null(comdist)) {
        comdist <- vegdist(comtable, method = vegdist.method)
    }

    # Auto-compute envdist from envmeta if not provided
    if (is.null(envdist)) {
        envdist <- dist(
            BBmisc::normalize(envmeta,
                method = "range",
                margin = 2
            ),
            method = "euclidean"
        )
    }

    # Function to call individual prediction methods
    call_method <- function(m) {
        result <- switch(m,
            "linearPred" = {
                linearPred(
                    comdist = comdist,
                    envdist = envdist,
                    sitenames = sitenames
                )
            },
            "mlPred" = {
                mlPred(
                    comdist = comdist,
                    envmeta = envmeta,
                    sitenames = sitenames
                )
            },
            "glmPred" = {
                glmPred(
                    comdist = comdist,
                    envmeta = envmeta,
                    sitenames = sitenames
                )
            },
            "gamPred" = {
                gamPred(
                    comtable = comtable,
                    envmeta = envmeta,
                    comdist = comdist,
                    sitenames = sitenames,
                    GAM.dist.method = GAM.dist.method
                )
            },
            "gdmPred" = {
                gdmPred(
                    comdist = comdist,
                    envmeta = envmeta,
                    sitenames = sitenames,
                    X = X,
                    Y = Y,
                    geo_enabled = geo_enabled
                )
            },
            "rfPred" = {
                rfPred(
                    comdist = comdist,
                    envmeta = envmeta,
                    sitenames = sitenames,
                    seed = seed
                )
            },
            "xgboostPred" = {
                xgboostPred(
                    comdist = comdist,
                    envmeta = envmeta,
                    sitenames = sitenames,
                    seed = seed,
                    params = xgboost.params
                )
            }
        )
        return(result)
    }

    # Execute method(s)
    if (length(method) == 1) {
        # Single method: return result vector directly
        return(call_method(method))
    } else {
        # Multiple methods: combine results into data frame
        results_list <- lapply(method, call_method)

        # Combine all results using cbind
        combined_result <- do.call(cbind, results_list)

        # Rename columns to indicate which method produced each result
        colnames(combined_result) <- method

        return(combined_result)
    }
}

#   #### testable data ####
#   # HSAUR3::birds, gardenflowers, watervoles
#
#   #### TODO: predict community then calculate diversity ####
#   # MicroEcoTools
#   # specificity
#   # microbiomeSeq


