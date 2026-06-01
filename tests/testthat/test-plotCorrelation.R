library(vegan)
data(varespec)
data(varechem)

test_that("plotCorrelation works with multiple methods result", {
  # Generate multiple methods result
  multi_result <- betaStability(
    comtable = varespec,
    envmeta = varechem,
    method = c("linearPred", "mlPred")
  )

  # Test that plotStability returns a ggplot object
  expect_s3_class(plotCorrelation(multi_result), "ggplot")
})
