test_that(
    "calcStability returns correct values",
    {
        # Test with known values
        expect_equal(calcStability(1, 0.5), 0.5)
        expect_equal(calcStability(0.5, 0.75), -0.5)
        expect_equal(calcStability(0.5, 0.5), 0)
        expect_equal(calcStability(0, 0), NaN)

        # Test with non-numeric inputs
        expect_error(calcStability("a", "b"))
        expect_error(calcStability(NA, NA))

        # Test with edge cases
        expect_equal(calcStability(1e-10, 1e-10), 0)
    }
)
