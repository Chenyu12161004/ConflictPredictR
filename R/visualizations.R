
#' Plot Feature Importance for the Random Forest Model
#'
#' Creates a bar plot of the most important features from the Random Forest model.
#'
#' @param top_n The number of top features to show. Defaults to 20.
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_minimal theme element_text
#' @importFrom dplyr as_tibble arrange desc top_n
#' @importFrom forcats fct_reorder
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @export
#' @examples
#' plot_rf_importance()
plot_rf_importance <- function(top_n = 20) {

  message("Loading Random Forest model...")
  model_rf <- readRDS(system.file("extdata", "final_rf_model.rds", package = "ConflictPredictR"))

  importance_df <- as.data.frame(model_rf$importance) %>%
    tibble::rownames_to_column("feature") %>%
    dplyr::arrange(dplyr::desc(MeanDecreaseGini))

  top_features <- dplyr::top_n(importance_df, top_n, MeanDecreaseGini)

  p_importance <- ggplot(top_features, aes(x = forcats::fct_reorder(feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(
      title = "Random Forest: Top Feature Importance",
      subtitle = "Based on Mean Decrease in Gini Impurity",
      x = "Feature",
      y = "Mean Decrease Gini"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold"))

  return(p_importance)
}

#' Plot Diagnostic Charts for the Random Forest Model
#'
#' Displays the Out-of-Bag (OOB) error rate curve and a confusion matrix of the model.
#'
#' The OOB error plot is static and reflects performance of the model during training.
#' The confusion matrix is dynamic: if newdata is provided, it shows performance
#' on that data; otherwise, it defaults to showing the performance on the OOB samples.
#'
#' @param newdata A data frame to evaluate for the confusion matrix. If NULL (the default),
#'   the function uses the internal Out-of-Bag confusion matrix of the model.
#' @return A combined ggplot object created by `grid.arrange`.
#' @importFrom ggplot2 ggplot aes geom_line labs theme_minimal theme geom_tile geom_text scale_fill_gradient scale_y_discrete element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom stats predict
#' @importFrom recipes bake
#' @importFrom dplyr select
#' @importFrom tidyr pivot_longer
#' @export
#' @examples
#' \dontrun{
#' # This shows the OOB error and the OOB confusion matrix
#' plot_rf_diagnostics()
#'
#' # This shows the OOB error and a new confusion matrix based on sample_variants
#' plot_rf_diagnostics(ConflictPredictR::sample_variants)
#' }
plot_rf_diagnostics <- function(newdata = NULL) {

  message("Loading Random Forest model...")
  model_rf <- readRDS(system.file("extdata", "final_rf_model.rds", package = "ConflictPredictR"))

  # OOB Error Plot
  if(!is.null(model_rf$err.rate)) {
    err_rate_df <- data.frame(
      trees = 1:model_rf$ntree,
      oob_error = model_rf$err.rate[, "OOB"],
      class0_error = model_rf$err.rate[, "0"],
      class1_error = model_rf$err.rate[, "1"]
    )

    p_error <- ggplot(err_rate_df, aes(x = trees)) +
      geom_line(aes(y = .data$oob_error, color = "OOB Error"), linewidth = 1) +
      geom_line(aes(y = .data$class0_error, color = "Class 0 Error"), linewidth = 0.7, alpha = 0.7) +
      geom_line(aes(y = .data$class1_error, color = "Class 1 Error"), linewidth = 0.7, alpha = 0.7) +
      labs(
        title = "Random Forest Error Rates (during training)",
        x = "Number of Trees", y = "Error Rate", color = "Error Type"
      ) +
      theme_minimal() + theme(legend.position = "bottom")
  } else {
    p_error <- ggplot() + labs(title = "OOB Error data not available.")
  }

  # Confusion Matrix Plot

  if (is.null(newdata)) {
    # No data provided. Default to showing the OOB confusion matrix.
    message("No 'newdata' provided. Displaying confusion matrix from OOB samples.")
    if (!is.null(model_rf$confusion)) {
      conf_matrix_df <- as.data.frame(model_rf$confusion)
      conf_matrix_df$Actual <- rownames(conf_matrix_df)
      conf_matrix_long <- tidyr::pivot_longer(conf_matrix_df, cols = c("0", "1"), names_to = "Predicted", values_to = "Freq")

      p_conf <- ggplot(conf_matrix_long, aes(x = Actual, y = Predicted, fill = Freq)) +
        geom_tile(color = "white") +
        geom_text(aes(label = Freq), color = "white", fontface = "bold", size = 6) +
        scale_fill_gradient(low = "lightblue", high = "darkblue") +
        labs(title = "Confusion Matrix (on OOB data)", x = "Actual Class", y = "Predicted Class", fill = "Count") +
        theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14))
    } else {
      p_conf <- ggplot() + labs(title = "OOB Confusion Matrix data not available.")
    }

  } else {
    # User provides newdata. Calculate a new confusion matrix.
    message("Calculating confusion matrix based on provided 'newdata'...")
    recipe_obj <- readRDS(system.file("extdata", "preprocessing_recipe.rds", package = "ConflictPredictR"))

    if (!"CLASS" %in% colnames(newdata)) {
      stop("The 'CLASS' column is missing from the input 'newdata'.")
    }
    actual_class <- newdata$CLASS

    baked_data <- recipes::bake(recipe_obj, new_data = newdata)
    pred_class <- stats::predict(model_rf, baked_data, type = "response")

    conf_matrix <- table(Predicted = pred_class, Actual = actual_class)
    conf_matrix_df <- as.data.frame(conf_matrix)

    p_conf <- ggplot(conf_matrix_df, aes(x = Actual, y = Predicted, fill = Freq)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Freq), color = "white", fontface = "bold", size = 6) +
      scale_fill_gradient(low = "lightblue", high = "darkblue") +
      labs(title = "Confusion Matrix (on provided data)", x = "Actual Class", y = "Predicted Class", fill = "Count") +
      theme_minimal() + theme(plot.title = element_text(face = "bold", size = 14))
  }

  # Combine and return the plots
  gridExtra::grid.arrange(p_error, p_conf, ncol = 2)
}
