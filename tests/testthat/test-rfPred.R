library(vegan)
data(varechem)
data(varespec)

test_that("random forest prediction no error", {
    expect_no_error(
        rfPred(
            vegdist(varespec, "bray"),
            varechem
        )
    )
})
