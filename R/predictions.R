
#' Predict Variant Pathogenicity using Random Forest
#'
#' This function takes new variant data and returns predictions from the
#' pre-trained Random Forest models.
#'
#' @param newdata A data frame of new variants to predict.
#' @return A tibble that includes the original data plus new columns for the predicted class and the probability of being pathogenic.
#' @importFrom recipes bake
#' @importFrom stats predict
#' @importFrom dplyr bind_cols select
#' @import randomForest
#' @export
#' @examples
#' # Run predictions on the built-in sample data
#' predictions <- predict_pathogenicity(PathoPrediction::sample_variants)
#' print(predictions)
predict_pathogenicity <- function(newdata) {

  if (nrow(newdata) == 0) {
    message("Input data has 0 rows. Returning an empty result.")
    return(
      dplyr::bind_cols(
        newdata,
        .pred_class = factor(character(0)),
        .pred_prob_pathogenic = numeric(0)
      )
    )
  }

  message("Loading Random Forest model and preprocessing recipe...")
  model_rf <- readRDS(system.file("extdata", "final_rf_model.rds", package = "PathoPrediction"))
  recipe_obj <- readRDS(system.file("extdata", "preprocessing_recipe.rds", package = "PathoPrediction"))

  message("Processing new data...")
  baked_data <- recipes::bake(recipe_obj, new_data = newdata)

  if ("CLASS" %in% colnames(baked_data)) {
    baked_data <- baked_data %>% dplyr::select(-"CLASS")
  }

  message("Making predictions...")
  # Predict probabilities
  pred_probs <- stats::predict(model_rf, baked_data, type = "prob")

  # Predict class
  pred_class <- stats::predict(model_rf, baked_data, type = "response")

  # Extract probability of being pathogenic (class "1")
  prob_pathogenic <- pred_probs[, "1"]

  # Combine results
  results <- dplyr::bind_cols(
    newdata,
    .pred_class = pred_class,
    .pred_prob_pathogenic = prob_pathogenic
  )

  message("Done.")
  return(results)
}
