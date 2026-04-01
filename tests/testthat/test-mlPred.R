data(varechem)
data(varespec)

test_that("multi linear prediction no error", {
  expect_no_error(
    mlPred(vegdist(varespec, "bray"),
                    varechem)
    )
})
