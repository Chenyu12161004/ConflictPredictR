# Test the validation function
test_that("validate_input_data works as expected", {
  expect_true(suppressWarnings(validate_input_data(sample_variants)))

  bad_data_missing_col <- sample_variants
  bad_data_missing_col$CADD_PHRED <- NULL
  expect_error(validate_input_data(bad_data_missing_col), "missing")
})

# Test the prediction function
test_that("predict_conflict_risk returns correct format", {
  risk_scores <- suppressWarnings(predict_conflict_risk(sample_variants))
  expect_s3_class(risk_scores, "data.frame")
  expect_true(all(c(".pred_conflict_class", ".pred_conflict_risk") %in% colnames(risk_scores)))
  expect_equal(nrow(risk_scores), nrow(sample_variants))

  empty_df <- sample_variants[0, ]
  empty_predictions <- predict_conflict_risk(empty_df)
  expect_equal(nrow(empty_predictions), 0)
})

# Test the plotting functions
test_that("plotting functions return ggplot objects", {
  p_importance <- plot_rf_importance()
  expect_s3_class(p_importance, "ggplot")

  p_diagnostics <- plot_rf_diagnostics()
  expect_s3_class(p_diagnostics, "gtable")
})
