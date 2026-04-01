test_that("beta diversity matrix same", {
  # Load example data
  data(varespec, package = "vegan")

  # Calculate beta diversity using your function
  actual <- betaDiv(varespec)

  # Calculate beta diversity using vegdist directly
  expected <- vegdist(varespec, "bray")

  # Strip the call attribute
  attr(actual, "call") <- NULL
  attr(expected, "call") <- NULL

  # Compare the results
  expect_equal(actual, expected)
})
