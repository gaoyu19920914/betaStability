#' Plot Stability Results
#'
#' This function takes the output of betaStability() and creates a point plot using ggplot2.
#'
#' @param stability_result The output from betaStability() function.
#' @param sitenames Optional vector of site names. If not provided, uses rownames from stability_result.
#'
#' @returns A ggplot2 plot object.
#'
#' @examples
#' data(varespec)
#' data(varechem)
#' results <- betaStability(comtable = varespec, envmeta = varechem, method = c("linearPred", "mlPred"))
#' plotStability(results)
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
plotStability <- function(stability_result, sitenames = NULL) {
  stability_result$site <- rownames(stability_result)
  # Convert to data frame if it's a vector
  # reshape to keep the method information
  df <- reshape2::melt(stability_result, id.vars = "site")
  colnames(df) <- c("site", "method", "stability")
  if (!is.null(sitenames)) {
    df$site <- sitenames
  }

  # Ensure site order is preserved
  df$site <- factor(df$site, levels = unique(df$site))

  # Create plot
  p <- ggplot(df, aes(x = site, y = stability, color = method)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ylim(-1, 1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Site", y = "Stability", color = "Method")

  return(p)
}
