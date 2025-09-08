# Test validation functions
library(testthat)
library(Seurat)

test_that("validateSeuratInput works correctly", {
  # Create mock Seurat object for testing
  skip_if_not_installed("Seurat")
  
  # Test with invalid input (not Seurat object)
  result <- validateSeuratInput("not_a_seurat_object")
  expect_false(result$valid)
  expect_true(length(result$errors) > 0)
  expect_match(result$errors[1], "not a Seurat object")
})

test_that("validateTLSParameters works correctly", {
  # Test valid parameters
  result <- validateTLSParameters(
    exp_threshold = 0.98,
    min_spots = 3,
    max_distance = 5,
    distance_multiplier = 2,
    alpha = 0.05
  )
  expect_true(result$valid)
  expect_equal(length(result$errors), 0)
  
  # Test invalid exp_threshold
  result <- validateTLSParameters(exp_threshold = 1.5)
  expect_false(result$valid)
  expect_true(any(grepl("exp_threshold must be between 0 and 1", result$errors)))
  
  # Test invalid min_spots
  result <- validateTLSParameters(min_spots = -1)
  expect_false(result$valid)
  expect_true(any(grepl("min_spots must be a positive integer", result$errors)))
  
  # Test invalid alpha
  result <- validateTLSParameters(alpha = 1.5)
  expect_false(result$valid)
  expect_true(any(grepl("alpha must be between 0 and 1", result$errors)))
  
  # Test warning for low exp_threshold
  result <- validateTLSParameters(exp_threshold = 0.7)
  expect_true(result$valid)  # Still valid, just warning
  expect_true(any(grepl("may result in many false positives", result$warnings)))
})