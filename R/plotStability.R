#' Plot Stability Results
#'
#' Creates a point plot of the stability results produced by
#' \code{betaStability()}. If a \code{coords} argument is supplied, a map of
#' the sites is drawn instead (optionally with an elevation raster as
#' background).
#'
#' @param stability_result The output from \code{betaStability()}.
#' @param coords Optional data frame of site coordinates. Columns may be named
#'   \code{x}/\code{y}, \code{X}/\code{Y}, \code{lon}/\code{lat}, or
#'   \code{longitude}/\code{latitude}; they are normalised internally. When
#'   supplied, \code{plotStability()} behaves like \code{plotStabilityMap()}.
#' @param sitenames Optional vector of site names. If not provided, uses
#'   rownames from \code{stability_result}. Users shall make sure the provided
#'   sitenames correspond to the rownames of \code{stability_result} in the
#'   correct order.
#' @param elev Logical. Only used when \code{coords} is supplied. If
#'   \code{TRUE}, an elevation raster is drawn behind the points using
#'   \pkg{elevatr}.
#'
#' @returns A \pkg{ggplot2} plot object.
#'
#' @examples
#' library(vegan)
#' library(ggplot2)
#' data(varespec)
#' data(varechem)
#' data(varecoords)
#' results <- betaStability(
#'     comtable = varespec,
#'     envmeta = varechem, method = c("linearPred", "mlPred")
#' )
#' plotStability(results)
#' plotStability(results, coords = varecoords)
#' \dontrun{
#' plotStability(results, coords = varecoords, elev = TRUE)
#' }
#'
#' # Alias
#' plotStabilityMap(results, varecoords)
#'
#' @importFrom reshape2 melt
#' @importFrom elevatr get_elev_raster
#' @importFrom grDevices terrain.colors
#' @importFrom raster as.data.frame
#' @import ggplot2
#' @export
plotStability <- function(stability_result,
                          coords = NULL,
                          sitenames = NULL,
                          elev = FALSE) {

  ## --- common prep ------------------------------------------------------
  stability_result$site <- rownames(stability_result)

  df <- reshape2::melt(stability_result, id.vars = "site")
  colnames(df) <- c("site", "method", "stability")

  if (!is.null(sitenames)) {
    df$site <- sitenames
  }

  ## --- case 1: no coordinates -> point plot -----------------------------
  if (is.null(coords)) {
    df$site <- factor(df$site, levels = unique(df$site))

    p <- ggplot(df, aes(x = site, y = stability, color = method)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      ylim(-1, 1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      labs(x = "Site", y = "Stability", color = "Method") +
      theme(legend.position = "bottom") +
      theme_bw()

    return(p)
  }

  ## --- case 2: coordinates supplied -> map ------------------------------
  ## Make sure coords carries a `site` column so merge() works, even if
  ## the user only supplies rownames.
  if (!"site" %in% colnames(coords)) {
    coords$site <- rownames(coords)
  }

  if (!all(stability_result$site %in% coords$site))
    stop("coordinates not complete or missing.")

  df <- merge(stability_result, coords, by = "site", all.x = TRUE)

  ## --- Normalise coordinate column names --------------------------------
  ## Accept x/y, X/Y (from sf::st_coordinates), lon/lat, longitude/latitude.
  nm <- names(df)
  x_candidates <- c("x", "X", "lon", "long", "longitude")
  y_candidates <- c("y", "Y", "lat", "latitude")

  x_col <- intersect(x_candidates, nm)
  y_col <- intersect(y_candidates, nm)

  if (length(x_col) == 0L || length(y_col) == 0L) {
    stop("`coords` must contain coordinate columns named ",
         "'x'/'y', 'X'/'Y', 'lon'/'lat', or 'longitude'/'latitude'.")
  }

  names(df)[names(df) == x_col[1L]] <- "x"
  names(df)[names(df) == y_col[1L]] <- "y"

  df <- reshape2::melt(df, id.vars = c("site", "x", "y"))
  colnames(df) <- c("site", "x", "y", "method", "stability")

  if (!is.null(sitenames)) {
    df$site <- sitenames
  }

  lim <- max(abs(df$stability), na.rm = TRUE)

  if (elev) {
    xlim_range <- range(df$x, na.rm = TRUE) + c(-1, 1)
    ylim_range <- range(df$y, na.rm = TRUE) + c(-1, 1)

    locations <- df[, c("x", "y")]
    elev_raster <- elevatr::get_elev_raster(
      locations = locations,
      z          = 6,
      clip       = "bbox",
      prj        = "EPSG:4326"
    )
    elev_df <- raster::as.data.frame(elev_raster, xy = TRUE, na.rm = TRUE)
    names(elev_df) <- c("x", "y", "elevation")

    p_map <- ggplot() +
      geom_raster(data = elev_df,
                  aes(x = x, y = y, fill = elevation),
                  alpha = 1) +
      scale_fill_gradientn(
        colors = terrain.colors(10),
        guide  = "none"
      ) +
      geom_point(data = df,
                 aes(x = x, y = y, color = stability),
                 size = 2, alpha = 1) +
      geom_text(data = df,
                aes(x = x, y = y, label = site),
                size = 2.5, vjust = -0.8,
                check_overlap = TRUE,
                show.legend = FALSE) +
      scale_color_gradient2(
        low      = "red",
        mid      = "white",
        high     = "blue",
        midpoint = 0,
        limits   = c(-lim, lim),
        name     = "betaStability",
        guide    = guide_colorbar(barwidth  = 1,
                                  barheight = 10,
                                  ticks     = TRUE)
      ) +
      facet_wrap(~ method) +
      coord_sf(xlim = xlim_range,
               ylim = ylim_range,
               expand = FALSE) +
      labs(x = "Longitude", y = "Latitude") +
      theme_bw() +
      theme(
        legend.position   = "right",
        strip.background  = element_rect(fill = "grey92", color = NA),
        strip.text        = element_text(face = "bold"),
        panel.grid.minor  = element_blank()
      )
    return(p_map)
  }

  p_scatter <- ggplot(df, aes(x = x, y = y, color = stability)) +
    geom_point(size = 2, alpha = 1) +
    geom_text(data = df,
              aes(x = x, y = y, label = site),
              size = 2.5, vjust = -0.8,
              check_overlap = TRUE,
              show.legend = FALSE) +
    scale_color_gradient2(
      low      = "red",
      mid      = "white",
      high     = "blue",
      midpoint = 0,
      limits   = c(-lim, lim),
      name     = "betaStability",
      guide    = guide_colorbar(barwidth  = 1,
                                barheight = 10,
                                ticks     = TRUE)
    ) +
    facet_wrap(~ method) +
    coord_fixed() +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      legend.position  = "right",
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  return(p_scatter)
}


#' @rdname plotStability
#' @export
plotStabilityMap <- function(stability_result,
                             coords,
                             sitenames = NULL,
                             elev = FALSE) {
  plotStability(stability_result,
                coords = coords,
                sitenames = sitenames,
                elev = elev)
}
