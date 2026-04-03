library(vegan)
data(varechem)
data(varespec)

test_that("generalized linear prediction no error", {
    expect_no_error(
        glmPred(
            vegdist(varespec, "bray"),
            varechem
        )
    )
})
