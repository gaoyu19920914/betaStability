library(vegan)
data(varechem)
data(varespec)
data(BCI)
data(BCI.env)

test_that("nonlinear prediction no error", {
    expect_no_error(
        gdmPred(
            vegdist(varespec, "bray"),
            varechem
        )
    )
})

test_that("gdmPred with X and Y coordinates", {
    expect_no_error(
        gdmPred(vegdist(BCI, "bray"),
            BCI.env[, c("Precipitation", "Elevation", "EnvHet")],
            X = BCI.env$UTM.EW,
            Y = BCI.env$UTM.NS
        )
    )
})
