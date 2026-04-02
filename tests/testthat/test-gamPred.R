library(testthat)
library(betaStability)
library(vegan)
data(varechem)
data(varespec)
library(mgcv)

test_that("generalized additive model prediction no error", {
    expect_no_error(
        gamPred(
            varespec,
            varechem
        )
    )
})
