library(testthat)
library(betaStability)
library(vegan)
data(varechem)
data(varespec)

test_that("nonlinear prediction no error", {
  expect_no_error(
    gdmPred(vegdist(varespec, "bray"),
            varechem)
  )
})
