library(testthat)

# Test data setup
data(varespec)
data(varechem)

# Test 1: Single method result
test_that("plotStability works with single method result", {
    # Generate single method result
    single_result <- betaStability(comtable = varespec, envmeta = varechem, method = "linearPred")

    # Test that plotStability returns a ggplot object
    expect_s3_class(plotStability(single_result), "ggplot")

    # Test with custom sitenames
    custom_sitenames <- paste("Site", 1:nrow(varespec))
    expect_s3_class(plotStability(single_result, sitenames = custom_sitenames), "ggplot")
})

# Test 2: Multiple methods result
test_that("plotStability works with multiple methods result", {
    # Generate multiple methods result
    multi_result <- betaStability(
        comtable = varespec, envmeta = varechem,
        method = c("linearPred", "mlPred")
    )

    # Test that plotStability returns a ggplot object
    expect_s3_class(plotStability(multi_result), "ggplot")

    # Test with custom sitenames
    custom_sitenames <- paste("Site", 1:nrow(varespec))
    expect_s3_class(plotStability(multi_result, sitenames = custom_sitenames), "ggplot")
})
