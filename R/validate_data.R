#' Validate Input Data for Prediction
#'
#' Checks if a new dataset has the structure and data types required by the
#' prediction model in `PathoPrediction`.
#'
#' @param newdata A data frame of new variants to validate.
#' @return Returns `TRUE` (invisibly) if the data is valid. If validation fails,
#'   it stops and prints an informative error message.
#' @export
#' @examples
#' # This will pass validation and print success messages
#' try(validate_input_data(PathoPrediction::sample_variants))
#'
#' # This example demonstrates a failure due to a missing column
#' bad_data <- PathoPrediction::sample_variants
#' bad_data$CADD_PHRED <- NULL # Remove a required column
#' try(validate_input_data(bad_data))
validate_input_data <- function(newdata) {

  message("Starting Data Validation")

  # Load the internal recipe of package to know the requirements
  recipe_obj <- readRDS(system.file("extdata", "preprocessing_recipe.rds", package = "PathoPrediction"))

  original_vars <- recipe_obj$var_info$variable

  # Is the input a data frame
  if (!is.data.frame(newdata)) {
    stop("Input must be a data.frame.")
  }

  # Are all required columns present
  missing_cols <- setdiff(original_vars, c("CLASS", colnames(newdata)))
  if (length(missing_cols) > 0) {
    stop(paste("The following required columns are missing:",
               paste(shQuote(missing_cols), collapse = ", ")))
  }
  message("All required columns are present.")

  # Is the data compatible with the recipe
  # Unexpected factor levels and other issues that `bake()` would fail on
  tryCatch({
    recipes::bake(recipe_obj, new_data = newdata)
    message("Data is compatible with the model's preprocessing steps.")
  }, error = function(e) {
    stop(paste("Data is not compatible with the model. This often means there are unexpected values (e.g., non-numeric data in a numeric column, or new factor levels).\n  Original error:", e$message))
  })

  message("Validation Successful")
  invisible(TRUE)
}
