library(testthat)
library(betaStability)
library(vegan)
data(varechem)
data(varespec)
library(mgcv, quietly = TRUE)

test_that("generalized additive model prediction no error", {
    expect_no_error(
        gamPred(
            varespec,
            varechem
        )
    )
})
