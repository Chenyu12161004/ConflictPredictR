# run_combined_models.R

# 0. Load required packages
pkgs <- c("tidyverse", "data.table", "tidymodels", "xgboost", "lightgbm", 
          "SHAPforxgboost", "pROC", "PRROC", "here", "ggplot2", "randomForest")
need <- setdiff(pkgs, rownames(installed.packages()))
if(length(need)) {
  install.packages(need, repos = "https://cloud.r-project.org", dependencies = TRUE)
}
invisible(lapply(pkgs, library, character.only = TRUE))

# 1. Set working directory and output paths
setwd("C:/Users/DELL/Desktop/study/BIO215/project")
OUT_XGB <- file.path("output/ml")
OUT_RF <- file.path("output/rf")
dir.create(OUT_XGB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_RF, recursive = TRUE, showWarnings = FALSE)
set.seed(215)

# 2. Read data
csv_path <- here("clinvar_conflicting.csv")
dt <- data.table::fread(csv_path)

# 3. Data cleaning - shared steps
dt <- dt %>% filter(CLASS %in% 0:1)

num_names <- names(dt)[sapply(dt, is.numeric)]
dt[, (num_names) := lapply(.SD, function(x) fifelse(is.na(x), median(x, na.rm = TRUE), x)),
   .SDcols = num_names]

chr_names <- names(dt)[sapply(dt, function(x) is.character(x) || is.factor(x))]
dt[, (chr_names) := lapply(.SD, function(x) {
  x <- factor(x)
  fct_explicit_na(x, na_level = "Missing")
}), .SDcols = chr_names]

# 4. Train/test split (stratified, 10% training)
data_split <- initial_split(dt, strata = CLASS, prop = 0.1)
train_df <- training(data_split) %>% mutate(CLASS = factor(CLASS, levels = c("0", "1")))
test_df  <- testing(data_split)  %>% mutate(CLASS = factor(CLASS, levels = c("0", "1")))

# 5. Design recipe (standardization, remove high correlation, zero variance)
rec <- recipe(CLASS ~ ., data = train_df) %>%
  step_rm(
    CHROM, REF, ALT, CLNDISDB, CLNDN, CLNHGVS, CLNVI, MC,
    Allele, Consequence, IMPACT, SYMBOL, Feature, EXON,
    cDNA_position, CDS_position, Protein_position, Amino_acids,
    SIFT, PolyPhen
  ) %>%
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  step_nzv(all_predictors()) %>%
  step_corr(all_numeric_predictors(), threshold = 0.95) %>%
  step_normalize(all_numeric_predictors())

# ==================== XGBoost Model Section ====================

# 6. 5-fold cross validation
cv_folds <- vfold_cv(train_df, v = 5, strata = CLASS)

# 7. XGBoost model definition
xgb_spec <-
  boost_tree(
    trees        = tune(),
    tree_depth   = tune(),
    learn_rate   = tune(),
    loss_reduction = tune(),
    min_n        = tune(),
    sample_size  = 0.8
  ) %>%
  set_engine("xgboost",
             objective   = "binary:logistic",
             eval_metric = "auc",
             early_stopping_rounds = 50,
             verbosity   = 0) %>%
  set_mode("classification")

# 8. workflow + grid tuning
xgb_wf <- workflow() %>% add_model(xgb_spec) %>% add_recipe(rec)

set.seed(215)
xgb_grid <- grid_space_filling(
  trees(), tree_depth(), learn_rate(), loss_reduction(), min_n(),
  size = 10
)

xgb_res <- tune_grid(
  xgb_wf,
  resamples = cv_folds,
  grid      = xgb_grid,
  metrics   = metric_set(roc_auc)
)

# 9. Select best parameters and refit
xgb_final <- finalize_workflow(xgb_wf, select_best(xgb_res, metric = "roc_auc")) %>%
  fit(train_df)

# 10. Test set performance
pred_prob <- predict(xgb_final, test_df, type = "prob")$.pred_1
lbl <- pull(test_df, CLASS)

roc_obj <- pROC::roc(lbl, pred_prob)
auc_val <- as.numeric(roc_obj$auc)
pr_obj  <- PRROC::pr.curve(scores.class0 = pred_prob[lbl == 1],
                           scores.class1 = pred_prob[lbl == 0], curve = TRUE)

writeLines(
  sprintf("Test AUROC: %.4f\nTest AUPRC: %.4f", auc_val, pr_obj$auc.integral),
  file.path(OUT_XGB, "performance.txt")
)

# 11. SHAP Interpretation
xgb_native <- extract_fit_parsnip(xgb_final)$fit

baked_data <- bake(prep(rec), train_df)
feat_names <- setdiff(colnames(baked_data), "CLASS")

batch_size <- 5000
n_batches <- ceiling(nrow(baked_data) / batch_size)
total_samples <- nrow(baked_data)

all_shap_scores <- matrix(0, nrow = total_samples, ncol = length(feat_names))
colnames(all_shap_scores) <- feat_names

df_for_shap <- baked_data %>%
  select(-CLASS) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.data.frame()

for (i in 1:n_batches) {
  start_idx <- (i-1) * batch_size + 1
  end_idx <- min(i * batch_size, total_samples)
  
  batch_data <- df_for_shap[start_idx:end_idx, , drop = FALSE]
  batch_matrix <- as.matrix(batch_data)
  
  shap_contrib <- predict(
    xgb_native,
    newdata = batch_matrix,
    predcontrib = TRUE,
    approxcontrib = FALSE,
    nthread = 1,
    strict_shape = FALSE
  )
  
  if (!is.null(shap_contrib)) {
    batch_shap <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]
    all_shap_scores[start_idx:end_idx, ] <- batch_shap
  }
  
  if (i %% 5 == 0) gc()
}

mean_shap_score <- colMeans(abs(all_shap_scores))

shap_values <- list(
  shap_score = all_shap_scores,
  BIAS0 = 0,
  mean_shap_score = mean_shap_score
)

xgb_mat <- bake(prep(rec), train_df) %>%
  select(-CLASS) %>%
  select(all_of(feat_names)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

if (!identical(colnames(shap_values$shap_score), colnames(xgb_mat))) {
  common_cols <- intersect(colnames(shap_values$shap_score), colnames(xgb_mat))
  
  if (length(common_cols) == 0) {
    shap_values$shap_score <- shap_values$shap_score[, 1:min(ncol(shap_values$shap_score), ncol(xgb_mat))]
    xgb_mat <- xgb_mat[, 1:min(ncol(shap_values$shap_score), ncol(xgb_mat))]
    
    colnames(shap_values$shap_score) <- paste0("Feature_", 1:ncol(shap_values$shap_score))
    colnames(xgb_mat) <- paste0("Feature_", 1:ncol(xgb_mat))
  } else {
    shap_values$shap_score <- shap_values$shap_score[, common_cols, drop = FALSE]
    xgb_mat <- xgb_mat[, common_cols, drop = FALSE]
  }
}

shap_long <- as.data.frame(shap_values$shap_score) %>%
  mutate(sample_id = 1:n()) %>%
  pivot_longer(
    cols = -sample_id,
    names_to = "feature",
    values_to = "shap_value"
  ) %>%
  left_join(
    as.data.frame(xgb_mat) %>% mutate(sample_id = 1:n()) %>%
      pivot_longer(
        cols = -sample_id,
        names_to = "feature", 
        values_to = "feature_value"
      ),
    by = c("sample_id", "feature")
  )

feature_order <- shap_long %>%
  group_by(feature) %>%
  summarise(importance = mean(abs(shap_value))) %>%
  arrange(importance) %>%
  pull(feature)

p_enhanced <- ggplot(shap_long, aes(x = factor(feature, levels = feature_order), y = shap_value)) +
  geom_point(aes(color = feature_value), alpha = 0.3, size = 1, position = position_jitter(width = 0.2)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  coord_flip() +
  scale_color_gradient2(
    low = "blue",  
    mid = "lightgreen",  
    high = "red", 
    midpoint = median(shap_long$feature_value, na.rm = TRUE),
    name = "Feature value" 
  )  + 
  labs(x = "Feature", 
       y = "SHAP Value", 
       title = "SHAP Feature Importance and Influence Direction",
       subtitle = "Color represents feature value magnitude, point position represents SHAP value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(OUT_XGB, "shap_summary.png"), p_enhanced, width = 12, height = 10, dpi = 300)

plot_fixed_waterfall <- function(shap_values, xgb_mat, row_id, top_n = 8) {
  
  baseline <- ifelse(is.list(shap_values$BIAS0), 
                     as.numeric(shap_values$BIAS0[[1]]), 
                     as.numeric(shap_values$BIAS0))
  
  current_shap <- as.data.frame(shap_values$shap_score[row_id, , drop = FALSE])
  current_feat <- as.data.frame(xgb_mat[row_id, , drop = FALSE])
  
  shap_vector <- as.numeric(current_shap[1, ])
  feat_vector <- as.numeric(current_feat[1, ])
  feature_names <- colnames(current_shap)
  
  final_value <- baseline + sum(shap_vector)
  
  feat_importance <- abs(shap_vector)
  top_indices <- order(feat_importance, decreasing = TRUE)[1:min(top_n, length(feature_names))]
  
  features <- feature_names[top_indices]
  shap_vals <- shap_vector[top_indices]
  feat_vals <- feat_vector[top_indices]
  
  plot_df <- data.frame(
    step = c("Baseline", features, "Final"),
    value = c(baseline, shap_vals, final_value),
    type = c("baseline", rep("feature", length(features)), "final"),
    label = c("Baseline", 
              paste0(features, "\n(", round(feat_vals, 2), ")"), 
              "Final Prediction"),
    stringsAsFactors = FALSE
  )
  
  plot_df$cumulative <- cumsum(plot_df$value)
  plot_df$start <- c(0, plot_df$cumulative[-nrow(plot_df)])
  plot_df$end <- plot_df$cumulative
  
  plot_df$color <- ifelse(plot_df$type == "baseline", "gray60",
                          ifelse(plot_df$type == "final", "black",
                                 ifelse(plot_df$value > 0, "red", "blue")))
  
  ggplot(plot_df, aes(x = factor(label, levels = label))) +
    geom_rect(aes(xmin = seq_along(label) - 0.4, 
                  xmax = seq_along(label) + 0.4,
                  ymin = start, ymax = end, fill = color),
              color = "black", linewidth = 0.3) +
    geom_text(aes(x = label, y = (start + end)/2, 
                  label = round(value, 3)), 
              size = 3, color = "white", fontface = "bold") +
    scale_fill_identity() +
    labs(
      title = paste("SHAP Waterfall Plot - Sample", row_id),
      subtitle = paste("Baseline:", round(baseline, 3), "Final:", round(final_value, 3)),
      x = "Features", 
      y = "Cumulative Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.title = element_text(face = "bold", size = 12)
    )
}

conf_rows <- which(train_df$CLASS == 1)[1:3]

for(i in seq_along(conf_rows)){
  tryCatch({
    current_row <- conf_rows[i]
    
    p_wp <- plot_fixed_waterfall(shap_values, xgb_mat, row_id = current_row)
    
    ggsave(file.path(OUT_XGB, sprintf("waterfall_conflict_%d.png", i)),
           p_wp, width = 12, height = 7, dpi = 300)
    
  }, error = function(e) {
  })
}

# 12. Create prediction function and export
predict_conflict_xgb <- function(newdata){
  newdata <- newdata %>% select(-CLASS)
  predict(xgb_final, newdata, type = "prob")$.pred_1
}

saveRDS(xgb_final, file.path(OUT_XGB, "final_xgb_model.rds"))
saveRDS(predict_conflict_xgb, file.path(OUT_XGB, "predict_conflict_fn.rds"))

# ==================== Random Forest Model Section ====================

# 13. Preprocess data (for randomForest)
prepped <- prep(rec, training = train_df)
train_processed <- bake(prepped, train_df) %>%
  mutate(CLASS = as.factor(CLASS))
test_processed <- bake(prepped, test_df) %>%
  mutate(CLASS = as.factor(CLASS))

# 14. Train Random Forest model (using empirical parameters)
n_features <- ncol(train_processed) - 1

best_params <- list(
  mtry = floor(sqrt(n_features)),
  ntree = 500,
  nodesize = 5
)

x_train <- train_processed %>% select(-CLASS)
y_train <- train_processed$CLASS
x_test <- test_processed %>% select(-CLASS)
y_test <- test_processed$CLASS

rf_final <- randomForest(
  x = x_train,
  y = y_train,
  ntree = best_params$ntree,
  mtry = best_params$mtry,
  nodesize = best_params$nodesize,
  importance = TRUE,
  do.trace = 50,
  proximity = FALSE,
  keep.forest = TRUE
)

# 15. Test set performance evaluation
pred_prob_rf <- predict(rf_final, x_test, type = "prob")[, "1"]

roc_obj_rf <- pROC::roc(y_test, pred_prob_rf)
auc_val_rf <- as.numeric(roc_obj_rf$auc)
pr_obj_rf  <- PRROC::pr.curve(scores.class0 = pred_prob_rf[y_test == "1"],
                              scores.class1 = pred_prob_rf[y_test == "0"], curve = TRUE)

writeLines(
  sprintf("Best Parameters: mtry=%d, ntree=%d, nodesize=%d\nTest AUROC: %.4f\nTest AUPRC: %.4f",
          best_params$mtry, best_params$ntree, best_params$nodesize, 
          auc_val_rf, pr_obj_rf$auc.integral),
  file.path(OUT_RF, "performance.txt")
)

# 16. Feature importance analysis
importance_df <- as.data.frame(rf_final$importance) %>%
  mutate(feature = rownames(.)) %>%
  arrange(desc(MeanDecreaseGini))

top_features <- head(importance_df, 20)

p_importance <- ggplot(top_features, aes(x = reorder(feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(
    title = "Random Forest Feature Importance (Top 20)",
    subtitle = "Based on Mean Decrease in Gini Impurity",
    x = "Feature",
    y = "Mean Decrease Gini"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray50")
  )

ggsave(file.path(OUT_RF, "feature_importance.png"), p_importance, width = 10, height = 8, dpi = 300)
write.csv(importance_df, file.path(OUT_RF, "feature_importance.csv"), row.names = FALSE)

# 17. Error curve (OOB error vs number of trees)
if(!is.null(rf_final$err.rate)) {
  err_rate_df <- data.frame(
    trees = 1:rf_final$ntree,
    oob_error = rf_final$err.rate[, "OOB"],
    class0_error = rf_final$err.rate[, "0"],
    class1_error = rf_final$err.rate[, "1"]
  )
  
  p_error <- ggplot(err_rate_df, aes(x = trees)) +
    geom_line(aes(y = oob_error, color = "OOB Error"), size = 1) +
    geom_line(aes(y = class0_error, color = "Class 0 Error"), size = 0.7, alpha = 0.7) +
    geom_line(aes(y = class1_error, color = "Class 1 Error"), size = 0.7, alpha = 0.7) +
    labs(
      title = "Random Forest Error Rates",
      x = "Number of Trees",
      y = "Error Rate",
      color = "Error Type"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(OUT_RF, "error_rates.png"), p_error, width = 10, height = 6, dpi = 300)
}

# 18. Confusion matrix
pred_class <- predict(rf_final, x_test)
conf_matrix <- table(Predicted = pred_class, Actual = y_test)

conf_matrix_df <- as.data.frame(conf_matrix)

p_conf <- ggplot(conf_matrix_df, aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "white", fontface = "bold", size = 6) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Confusion Matrix",
    x = "Actual Class",
    y = "Predicted Class",
    fill = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 10)
  )

ggsave(file.path(OUT_RF, "confusion_matrix.png"), p_conf, width = 8, height = 6, dpi = 300)

accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
sensitivity <- conf_matrix[2, 2] / sum(conf_matrix[, 2])
specificity <- conf_matrix[1, 1] / sum(conf_matrix[, 1])
precision <- conf_matrix[2, 2] / sum(conf_matrix[2, ])

metrics_df <- data.frame(
  Metric = c("Accuracy", "Sensitivity (Recall)", "Specificity", "Precision", "AUROC", "AUPRC"),
  Value = c(accuracy, sensitivity, specificity, precision, auc_val_rf, pr_obj_rf$auc.integral)
)

write.csv(metrics_df, file.path(OUT_RF, "evaluation_metrics.csv"), row.names = FALSE)

# 19. Create prediction function
predict_conflict_rf <- function(newdata){
  newdata_processed <- bake(prepped, newdata)
  
  if("CLASS" %in% colnames(newdata_processed)) {
    newdata_processed <- newdata_processed %>% select(-CLASS)
  }
  
  pred_probs <- predict(rf_final, newdata_processed, type = "prob")
  
  if(ncol(pred_probs) == 2 && "1" %in% colnames(pred_probs)) {
    return(pred_probs[, "1"])
  } else if(ncol(pred_probs) == 2) {
    return(pred_probs[, 2])
  } else {
    return(pred_probs)
  }
}

# 20. Save model and preprocessing objects
saveRDS(rf_final, file.path(OUT_RF, "final_rf_model.rds"))
saveRDS(prepped, file.path(OUT_RF, "preprocessing_recipe.rds"))
saveRDS(predict_conflict_rf, file.path(OUT_RF, "predict_conflict_rf_fn.rds"))