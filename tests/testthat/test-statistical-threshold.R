# Test statistical threshold calculation
library(testthat)

test_that("calculateStatisticalThreshold works correctly", {
  # Create test data
  set.seed(123)
  expression_values <- c(rnorm(100, mean = 0, sd = 1), rnorm(20, mean = 3, sd = 0.5))
  
  # Test basic functionality
  result <- calculateStatisticalThreshold(expression_values, alpha = 0.05, method = "fdr")
  
  expect_type(result, "list")
  expect_true(all(c("threshold", "method", "p_values", "adjusted_p_values", "n_significant", "alpha") %in% names(result)))
  expect_true(is.numeric(result$threshold))
  expect_true(result$threshold > 0)
  expect_equal(result$alpha, 0.05)
  expect_match(result$method, "statistical_fdr")
  
  # Test with small sample size
  small_sample <- c(1, 2, 3, 4, 5)
  result_small <- calculateStatisticalThreshold(small_sample)
  
  expect_true(result_small$method == "quantile_fallback")
  expect_true(all(is.na(result_small$p_values)))
  
  # Test with different methods
  result_bonferroni <- calculateStatisticalThreshold(expression_values, method = "bonferroni")
  expect_match(result_bonferroni$method, "statistical_bonferroni")
  
  # Bonferroni should be more conservative (fewer significant spots)
  expect_lte(result_bonferroni$n_significant, result$n_significant)
  
  # Test with missing values
  expression_with_na <- c(expression_values, NA, NA, NA)
  result_na <- calculateStatisticalThreshold(expression_with_na)
  expect_type(result_na, "list")
  expect_true(is.finite(result_na$threshold))
  
  # Test edge cases
  all_same <- rep(5, 100)
  result_same <- calculateStatisticalThreshold(all_same)
  expect_equal(result_same$n_significant, 0)  # No variation, no significant spots
})

test_that("calculateStatisticalThreshold parameter validation", {
  test_data <- rnorm(100)
  
  # Test invalid alpha values
  expect_error(calculateStatisticalThreshold(test_data, alpha = -0.1), NA)  # Should not error, but may give warning
  expect_error(calculateStatisticalThreshold(test_data, alpha = 1.1), NA)   # Should not error, but may give warning
  
  # Test invalid method
  expect_error(calculateStatisticalThreshold(test_data, method = "invalid_method"))
})