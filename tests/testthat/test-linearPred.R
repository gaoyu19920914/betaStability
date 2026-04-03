library(vegan)
library(stats)
data(varechem)
data(varespec)

test_that("linear prediction no error", {
    expect_no_error(
        linearPred(
            vegdist(varespec, "bray"),
            dist(
                BBmisc::normalize(varechem,
                    method = "range",
                    margin = 2
                ),
                method = "euclidean"
            )
        )
    )
})
