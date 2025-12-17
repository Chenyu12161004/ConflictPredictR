# run_optimized_rf.R

# 0. Load required packages
pkgs <- c("tidyverse", "data.table", "tidymodels", "recipes", "pROC", 
          "PRROC", "here", "ggplot2", "randomForest", "caret")
need <- setdiff(pkgs, rownames(installed.packages()))
if(length(need)) {
  install.packages(need, repos = "https://cloud.r-project.org", dependencies = TRUE)
}
invisible(lapply(pkgs, library, character.only = TRUE))

# 1. Set working directory and output path
setwd("C:/Users/DELL/Desktop/study/BIO215/project")
OUT_RF <- file.path("output/rf")
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

# 4. Train/test split (stratified, 80% training)
data_split <- initial_split(dt, strata = CLASS, prop = 0.8)
train_df <- training(data_split) %>% mutate(CLASS = factor(CLASS, levels = c("0", "1")))
test_df  <- testing(data_split)  %>% mutate(CLASS = factor(CLASS, levels = c("0", "1")))

# 5. Design optimized recipe
rec <- recipe(CLASS ~ ., data = train_df) %>%
  step_rm(
    CHROM, REF, ALT, CLNDISDB, CLNDN, CLNHGVS, CLNVI, MC,
    Allele, Consequence, IMPACT, SYMBOL, Feature, EXON,
    cDNA_position, CDS_position, Protein_position, Amino_acids,
    SIFT, PolyPhen
  ) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_unknown(all_nominal_predictors(), new_level = "Missing") %>%
  step_dummy(all_nominal_predictors(), one_hot = FALSE) %>%
  step_nzv(all_predictors()) %>%
  step_corr(all_numeric_predictors(), threshold = 0.9) %>%
  step_normalize(all_numeric_predictors())

# 6. Preprocess data
prepped <- prep(rec, training = train_df)
train_processed <- bake(prepped, train_df) %>%
  mutate(CLASS = as.factor(CLASS))
test_processed <- bake(prepped, test_df) %>%
  mutate(CLASS = as.factor(CLASS))

# 7. Optimize random forest parameters
x_train <- train_processed %>% select(-CLASS)
y_train <- train_processed$CLASS
x_test <- test_processed %>% select(-CLASS)
y_test <- test_processed$CLASS

n_features <- ncol(train_processed) - 1

param_grid <- expand.grid(
  mtry = c(floor(sqrt(n_features)), 
           floor(n_features/3), 
           floor(n_features/2)),
  ntree = c(500, 800, 1000),
  nodesize = c(1, 3, 5)
)

best_auc <- 0
best_params <- list()

for(i in 1:nrow(param_grid)) {
  set.seed(215 + i)
  folds <- createFolds(as.numeric(y_train), k = 5, list = TRUE, returnTrain = FALSE)
  
  fold_aucs <- numeric(5)
  
  for(fold in 1:5) {
    val_idx <- folds[[fold]]
    train_fold_x <- x_train[-val_idx, ]
    train_fold_y <- y_train[-val_idx]
    val_fold_x <- x_train[val_idx, ]
    val_fold_y <- y_train[val_idx]
    
    rf_fold <- randomForest(
      x = train_fold_x,
      y = train_fold_y,
      ntree = param_grid$ntree[i],
      mtry = param_grid$mtry[i],
      nodesize = param_grid$nodesize[i],
      importance = FALSE,
      do.trace = FALSE
    )
    
    pred_prob <- predict(rf_fold, val_fold_x, type = "prob")[, "1"]
    fold_aucs[fold] <- pROC::roc(val_fold_y, pred_prob)$auc
  }
  
  mean_auc <- mean(fold_aucs)
  
  if(mean_auc > best_auc) {
    best_auc <- mean_auc
    best_params <- list(
      mtry = param_grid$mtry[i],
      ntree = param_grid$ntree[i],
      nodesize = param_grid$nodesize[i]
    )
  }
}

# 8. Train final random forest model (with class weights)
classwt <- c("0" = 1, "1" = 2)

rf_final <- randomForest(
  x = x_train,
  y = y_train,
  ntree = best_params$ntree,
  mtry = best_params$mtry,
  nodesize = best_params$nodesize,
  classwt = classwt,
  importance = TRUE,
  do.trace = 50,
  proximity = FALSE,
  keep.forest = TRUE,
  strata = y_train
)

# 9. Test set performance evaluation
pred_prob_rf <- predict(rf_final, x_test, type = "prob")[, "1"]
roc_obj_rf <- pROC::roc(y_test, pred_prob_rf)
auc_val_rf <- as.numeric(roc_obj_rf$auc)
pr_obj_rf  <- PRROC::pr.curve(scores.class0 = pred_prob_rf[y_test == "1"],
                              scores.class1 = pred_prob_rf[y_test == "0"], curve = TRUE)

performance_text <- sprintf(
  "Best parameters: mtry=%d, ntree=%d, nodesize=%d\nClass weights: 0=%.1f, 1=%.1f\nTest AUROC: %.4f\nTest AUPRC: %.4f",
  best_params$mtry, best_params$ntree, best_params$nodesize,
  classwt["0"], classwt["1"],
  auc_val_rf, pr_obj_rf$auc.integral
)

writeLines(performance_text, file.path(OUT_RF, "performance.txt"))

# 10. Feature importance analysis
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

# 11. Error curve (OOB error vs number of trees)
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

# 12. Confusion matrix
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

# 13. ROC and PR curves
roc_data <- data.frame(
  fpr = 1 - roc_obj_rf$specificities,
  tpr = roc_obj_rf$sensitivities
)

p_roc <- ggplot(roc_data, aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  annotate("text", x = 0.7, y = 0.3, 
           label = sprintf("AUROC = %.4f", auc_val_rf),
           size = 5, color = "red") +
  labs(
    title = "ROC Curve",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_minimal()

ggsave(file.path(OUT_RF, "roc_curve.png"), p_roc, width = 8, height = 6, dpi = 300)

pr_data <- data.frame(
  recall = pr_obj_rf$curve[, 1],
  precision = pr_obj_rf$curve[, 2]
)

p_pr <- ggplot(pr_data, aes(x = recall, y = precision)) +
  geom_line(color = "red", size = 1) +
  annotate("text", x = 0.7, y = 0.3, 
           label = sprintf("AUPRC = %.4f", pr_obj_rf$auc.integral),
           size = 5, color = "blue") +
  labs(
    title = "Precision-Recall Curve",
    x = "Recall (Sensitivity)",
    y = "Precision"
  ) +
  theme_minimal()

ggsave(file.path(OUT_RF, "pr_curve.png"), p_pr, width = 8, height = 6, dpi = 300)

# 14. Prediction function wrapper
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

# 15. Save model and preprocessing objects
saveRDS(rf_final, file.path(OUT_RF, "final_rf_model.rds"))
saveRDS(prepped, file.path(OUT_RF, "preprocessing_recipe.rds"))
saveRDS(predict_conflict_rf, file.path(OUT_RF, "predict_conflict_rf_fn.rds"))