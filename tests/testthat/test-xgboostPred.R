data(varechem)
data(varespec)

test_that("xgboost prediction no error", {
    expect_no_error(
        xgboostPred(
            vegdist(varespec, "bray"),
            varechem
        )
    )
})
