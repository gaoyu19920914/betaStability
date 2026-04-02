data(varechem)
data(varespec)

test_that("betaStability single method (linearPred) works correctly", {
  result <- betaStability(comtable = varespec,
                          envmeta = varechem,
                          method = "linearPred")

  # Check no error
  expect_no_error(result)

  # Check result is a data frame
  expect_s3_class(result, "data.frame")

  # Check has exactly 1 column
  expect_equal(ncol(result), 1)

  # Check column name
  expect_equal(colnames(result), "stability_Linear")

  # Check number of rows matches number of sites
  expect_equal(nrow(result), nrow(varespec))
})

test_that("betaStability multiple methods works correctly", {
  result <- betaStability(comtable = varespec,
                          envmeta = varechem,
                          method = c("linearPred", "mlPred"))

  # Check no error
  expect_no_error(result)

  # Check result is a data frame
  expect_s3_class(result, "data.frame")

  # Check has exactly 2 columns
  expect_equal(ncol(result), 2)

  # Check column names match the methods
  expect_equal(colnames(result), c("linearPred", "mlPred"))

  # Check number of rows matches number of sites
  expect_equal(nrow(result), nrow(varespec))
})
