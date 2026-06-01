#' Plot Correlation Matrix
#'
#' This function takes a result dataframe from `betaStability()` and creates a
#' faceted scatter plot matrix to visualize correlations between different
#' stability quantification methods.
#'
#' @param data A dataframe containing stability results from `betaStability()`.
#'   Must have at least 2 numeric columns.
#' @param method Correlation method to use. Default is "spearman". Other
#'   options include "pearson" and "kendall".
#' @importFrom stats cor.test
#' @returns A ggplot2 plot object showing pairwise correlations between columns.
#'
#' @examples
#' library(vegan)
#' library(ggplot2)
#' data(varespec)
#' data(varechem)
#' results <- betaStability(
#'     comtable = varespec,
#'     envmeta = varechem,
#'     method = c("linearPred", "mlPred", "glmPred")
#' )
#' plotCorrelation(results)
#'
#' @import ggplot2
#' @export
plotCorrelation <- function(data, method = "spearman") {
  if (!is.data.frame(data) || ncol(data) < 2) {
    stop("data must be a dataframe with at least 2 columns")
  }

  numeric_cols <- names(data)[vapply(data, is.numeric, FUN.VALUE = logical(1))]

  if (length(numeric_cols) < 2) {
    stop("data must contain at least 2 numeric columns")
  }

  pair_data <- do.call(rbind, lapply(numeric_cols, function(yvar) {
    do.call(rbind, lapply(numeric_cols, function(xvar) {
      data.frame(
        xval = data[[xvar]],
        yval = data[[yvar]],
        xvar = xvar,
        yvar = yvar,
        stringsAsFactors = FALSE
      )
    }))
  }))

  cor_results <- do.call(rbind, lapply(numeric_cols, function(yvar) {
    do.call(rbind, lapply(numeric_cols, function(xvar) {
      cor_test <- cor.test(data[[xvar]], data[[yvar]], method = method)
      cor_label <- ifelse(method == "spearman", "rho", "r")
      data.frame(
        xvar = xvar,
        yvar = yvar,
        label = sprintf("%s = %.2f, p = %.3f",
                        cor_label,
                        cor_test$estimate,
                        cor_test$p.value),
        stringsAsFactors = FALSE
      )
    }))
  }))

  p <- ggplot(pair_data, aes(x = xval, y = yval)) +
    geom_point(color = "black", size = 1.5, alpha = 0.6) +
    geom_smooth(method = "loess", se = TRUE, formula = y ~ x) +
    geom_text(data = cor_results, aes(x = -Inf, y = Inf, label = label),
              hjust = 0, vjust = 1, size = 2.5, inherit.aes = FALSE) +
    facet_grid(rows = vars(yvar), cols = vars(xvar)) +
    labs(title = paste(toupper(method),
                       "Correlation of Stability Predictions between Methods"),
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(strip.text = element_text(size = 8),
          axis.text = element_text(size = 6))

  return(p)
}
