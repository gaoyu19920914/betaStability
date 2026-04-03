library(vegan)
library(mgcv, quietly = TRUE)
data(varechem)
data(varespec)

test_that("generalized additive model prediction no error", {
    expect_no_error(
        gamPred(
            varespec,
            varechem
        )
    )
})
