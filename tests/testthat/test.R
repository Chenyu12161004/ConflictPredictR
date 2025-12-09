
# Test the validation function
test_that("validate_input_data works as expected", {

  # A good dataset should pass
  expect_true(validate_input_data(sample_variants))

  # A dataset missing a column should fail
  bad_data_missing_col <- sample_variants
  bad_data_missing_col$CADD_PHRED <- NULL
  expect_error(validate_input_data(bad_data_missing_col), "missing")

  # A dataset with wrong data type should fail
  bad_data_wrong_type <- sample_variants
  bad_data_wrong_type$CADD_PHRED <- as.character(bad_data_wrong_type$CADD_PHRED)
  expect_error(validate_input_data(bad_data_wrong_type), "compatible")
})

# Test the prediction function
test_that("predict_pathogenicity returns correct format", {

  # Prediction on valid data should return a data frame with specific columns
  predictions <- predict_pathogenicity(sample_variants)
  expect_s3_class(predictions, "data.frame")
  expect_true(all(c(".pred_class", ".pred_prob_pathogenic") %in% colnames(predictions)))
  expect_equal(nrow(predictions), nrow(sample_variants))

  # Prediction on an empty data frame should return an empty data frame
  empty_df <- sample_variants[0, ]
  empty_predictions <- predict_pathogenicity(empty_df)
  expect_equal(nrow(empty_predictions), 0)
})

# Test the plotting functions (simple tests)
test_that("plotting functions return ggplot objects", {

  # Check that the functions run without error and return a ggplot object
  p_importance <- plot_rf_importance()
  expect_s3_class(p_importance, "ggplot")

  # The diagnostics plot returns a 'gtable' object from grid.arrange
  p_diagnostics <- plot_rf_diagnostics()
  expect_s3_class(p_diagnostics, "gtable")
})
