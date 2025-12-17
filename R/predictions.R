
#' Predict the Risk of Conflicting Interpretations for Variants
#'
#' This function takes new variant data and returns a "Conflict Risk Score"
#' from the pre-trained Random Forest model.
#'
#' @param newdata A data frame of new variants to predict.
#' @return A tibble that includes the original data plus new columns for the
#'   predicted conflict class (`.pred_conflict_class`) and the conflict risk
#'   score (`.pred_conflict_risk`).
#' @importFrom recipes bake
#' @importFrom stats predict
#' @importFrom dplyr bind_cols select
#' @import randomForest
#' @export
#' @examples
#' # Run predictions on the built-in sample data
#' risk_scores <- predict_conflict_risk(ConflictPredictR::sample_variants)
#' # View a subset of the results
#' print(head(risk_scores[, c("CADD_PHRED", ".pred_conflict_class", ".pred_conflict_risk")]))
predict_conflict_risk <- function(newdata) {

  if (nrow(newdata) == 0) {
    message("Input data has 0 rows. Returning an empty result.")
    return(
      dplyr::bind_cols(
        newdata,
        .pred_conflict_class = factor(character(0), levels = c("0", "1")),
        .pred_conflict_risk = numeric(0)
      )
    )
  }

  message("Loading conflict risk model and preprocessing recipe...")
  model_rf <- readRDS(system.file("extdata", "final_rf_model.rds", package = "ConflictPredictR"))
  recipe_obj <- readRDS(system.file("extdata", "preprocessing_recipe.rds", package = "ConflictPredictR"))

  message("Processing new data...")
  baked_data <- recipes::bake(recipe_obj, new_data = newdata)

  if ("CLASS" %in% colnames(baked_data)) {
    baked_data <- baked_data %>% dplyr::select(-"CLASS")
  }

  message("Calculating conflict risk scores...")
  pred_probs <- predict(model_rf, baked_data, type = "prob")
  pred_class <- predict(model_rf, baked_data, type = "response")
  prob_conflict <- pred_probs[, "1"]

  results <- dplyr::bind_cols(
    newdata,
    .pred_conflict_class = pred_class,
    .pred_conflict_risk = prob_conflict
  )

  message("Done.")
  return(results)
}
