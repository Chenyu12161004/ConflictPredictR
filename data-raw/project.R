# run_combined_models.R

# 0. 加载所有需要的包
pkgs <- c("tidyverse", "data.table", "tidymodels", "xgboost", "lightgbm", 
          "SHAPforxgboost", "pROC", "PRROC", "here", "ggplot2", "randomForest")
need <- setdiff(pkgs, rownames(installed.packages()))
if(length(need)) {
  cat("Installing missing packages:", paste(need, collapse = ", "), "\n")
  install.packages(need, repos = "https://cloud.r-project.org", dependencies = TRUE)
}
invisible(lapply(pkgs, library, character.only = TRUE))

# 1. 设置工作目录与输出路径
setwd("C:/Users/DELL/Desktop/study/BIO215/project")
OUT_XGB <- file.path("output/ml")
OUT_RF <- file.path("output/rf")
dir.create(OUT_XGB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_RF, recursive = TRUE, showWarnings = FALSE)
set.seed(215)  # 固定随机种子

# 2. 读取数据
csv_path <- here("clinvar_conflicting.csv")
dt <- data.table::fread(csv_path)

# 3. 数据清洗 - 共用部分
# 3.1 只保留 CLASS 为 0/1 的行
dt <- dt %>% filter(CLASS %in% 0:1)

# 3.2 数值列：缺失 -> 中位数
num_names <- names(dt)[sapply(dt, is.numeric)]
dt[, (num_names) := lapply(.SD, function(x) fifelse(is.na(x), median(x, na.rm = TRUE), x)),
   .SDcols = num_names]

# 3.3 字符/因子列：缺失 -> "Missing"
chr_names <- names(dt)[sapply(dt, function(x) is.character(x) || is.factor(x))]
dt[, (chr_names) := lapply(.SD, function(x) {
  x <- factor(x)
  fct_explicit_na(x, na_level = "Missing")
}), .SDcols = chr_names]

# 4. 训练/测试划分（分层抽样，10% 训练）
data_split <- initial_split(dt, strata = CLASS, prop = 0.1)
train_df <- training(data_split) %>% mutate(CLASS = factor(CLASS, levels = c("0", "1")))
test_df  <- testing(data_split)  %>% mutate(CLASS = factor(CLASS, levels = c("0", "1")))

# 5. 设计 recipe（标准化、去除高相关、零方差）
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

# ==================== XGBoost 模型部分 ====================
cat("\n=== Training XGBoost Model ===\n")

# 6. 5 折交叉验证
cv_folds <- vfold_cv(train_df, v = 5, strata = CLASS)

# 7. XGBoost 模型定义
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
             early_stop_rounds = 50,
             verbosity   = 0) %>%
  set_mode("classification")

# 8. workflow + 网格调参
xgb_wf <- workflow() %>% add_model(xgb_spec) %>% add_recipe(rec)

set.seed(215)
xgb_grid <- grid_latin_hypercube(
  trees(), tree_depth(), learn_rate(), loss_reduction(), min_n(),
  size = 10
)

xgb_res <- tune_grid(
  xgb_wf,
  resamples = cv_folds,
  grid      = xgb_grid,
  metrics   = metric_set(roc_auc)
)

# 9. 取最优重新拟合
xgb_final <- finalize_workflow(xgb_wf, select_best(xgb_res, metric = "roc_auc")) %>%
  fit(train_df)

# 10. 测试集性能
pred_prob <- predict(xgb_final, test_df, type = "prob")$.pred_1
lbl <- pull(test_df, CLASS)

roc_obj <- pROC::roc(lbl, pred_prob)
auc_val <- as.numeric(roc_obj$auc)
pr_obj  <- PRROC::pr.curve(scores.class0 = pred_prob[lbl == 1],
                           scores.class1 = pred_prob[lbl == 0], curve = TRUE)

# 保存指标
writeLines(
  sprintf("Test AUROC: %.4f\nTest AUPRC: %.4f", auc_val, pr_obj$auc.integral),
  file.path(OUT_XGB, "performance.txt")
)

cat(sprintf("XGBoost Test Set Performance:\nAUROC: %.4f\nAUPRC: %.4f\n", auc_val, pr_obj$auc.integral))

# 11. SHAP 解释
# 11.1 提取原生 xgb.Booster
xgb_native <- extract_fit_parsnip(xgb_final)$fit

# 11.2.1. 取出训练时真正用的特征顺序
feat_names <- extract_fit_parsnip(xgb_final)$fit$feature_names

# 11.2.2. 按同样顺序、同样类型转矩阵
xgb_mat <- bake(prep(rec), train_df) %>% 
  select(-CLASS) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  as.matrix() %>% 
  .[, feat_names, drop = FALSE]

# 11.2.3. 计算 SHAP
shap_values <- SHAPforxgboost::shap.values(
  xgb_model = extract_fit_parsnip(xgb_final)$fit,
  X_train   = xgb_mat
)

# 11.3 全局 summary 图
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

# 计算特征重要性排序
feature_order <- shap_long %>%
  group_by(feature) %>%
  summarise(importance = mean(abs(shap_value))) %>%
  arrange(importance) %>%
  pull(feature)

# 绘制增强版summary图
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
  labs(x = "feature ", 
       y = "SHAP value ", 
       title =" Importance and Influence Direction of SHAP Feature ", 
       subtitle = "Color represents the size of feature value, point position represents SHAP value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )


ggsave(file.path(OUT_XGB, "shap_summary.png"), p_enhanced, width = 12, height = 10, dpi = 300)

# 11.4 对前 3 个冲突样本画 waterfall
plot_fixed_waterfall <- function(shap_values, xgb_mat, row_id, top_n = 8) {
  
  # 获取基准值 - 处理可能是列表的情况
  baseline <- ifelse(is.list(shap_values$BIAS0), 
                     as.numeric(shap_values$BIAS0[[1]]), 
                     as.numeric(shap_values$BIAS0))
  
  # 获取当前样本数据
  current_shap <- as.data.frame(shap_values$shap_score[row_id, , drop = FALSE])
  current_feat <- as.data.frame(xgb_mat[row_id, , drop = FALSE])
  
  # 转换为数值向量
  shap_vector <- as.numeric(current_shap[1, ])
  feat_vector <- as.numeric(current_feat[1, ])
  feature_names <- colnames(current_shap)
  
  final_value <- baseline + sum(shap_vector)
  
  # 选择最重要的特征
  feat_importance <- abs(shap_vector)
  top_indices <- order(feat_importance, decreasing = TRUE)[1:min(top_n, length(feature_names))]
  
  # 构建数据
  features <- feature_names[top_indices]
  shap_vals <- shap_vector[top_indices]
  feat_vals <- feat_vector[top_indices]
  
  # 创建数据框
  plot_df <- data.frame(
    step = c("Baseline", features, "Final"),
    value = c(baseline, shap_vals, final_value),
    type = c("baseline", rep("feature", length(features)), "final"),
    label = c("Baseline", 
              paste0(features, "\n(", round(feat_vals, 2), ")"), 
              "Final Prediction"),
    stringsAsFactors = FALSE
  )
  
  # 计算累积位置
  plot_df$cumulative <- cumsum(plot_df$value)
  plot_df$start <- c(0, plot_df$cumulative[-nrow(plot_df)])
  plot_df$end <- plot_df$cumulative
  
  # 设置颜色
  plot_df$color <- ifelse(plot_df$type == "baseline", "gray60",
                          ifelse(plot_df$type == "final", "black",
                                 ifelse(plot_df$value > 0, "red", "blue")))
  
  # 绘图
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

# 先检查数据结构
cat("=== XGBoost SHAP Data Structure Check ===\n")
cat("SHAP BIAS0 type:", class(shap_values$BIAS0), "\n")
cat("SHAP BIAS0 length:", length(shap_values$BIAS0), "\n")
cat("SHAP BIAS0 value:", if(is.list(shap_values$BIAS0)) shap_values$BIAS0[[1]] else shap_values$BIAS0, "\n\n")

cat("SHAP scores dimensions:", dim(shap_values$shap_score), "\n")
cat("XGB matrix dimensions:", dim(xgb_mat), "\n")
cat("Train df CLASS sum:", sum(train_df$CLASS == 1), "\n\n")

# 使用修复版本
conf_rows <- which(train_df$CLASS == 1)[1:3]
cat("Processing XGBoost SHAP samples:", conf_rows, "\n")

for(i in seq_along(conf_rows)){
  tryCatch({
    current_row <- conf_rows[i]
    cat("Creating XGBoost waterfall plot for sample", current_row, "...\n")
    
    p_wp <- plot_fixed_waterfall(shap_values, xgb_mat, row_id = current_row)
    print(p_wp)
    
    ggsave(file.path(OUT_XGB, sprintf("waterfall_conflict_%d.png", i)),
           p_wp, width = 12, height = 7, dpi = 300)
    
    cat("Successfully created XGBoost waterfall plot for sample", current_row, "\n")
    
  }, error = function(e) {
    cat("Error with XGBoost sample", conf_rows[i], ":", e$message, "\n")
  })
}

# 12. 封装"一键预测函数"并导出
predict_conflict_xgb <- function(newdata){
  # newdata 必须包含与训练时相同的特征列（不含 CLASS）
  newdata <- newdata %>% select(-CLASS)
  predict(xgb_final, newdata, type = "prob")$.pred_1
}

# 保存模型 & 函数
saveRDS(xgb_final, file.path(OUT_XGB, "final_xgb_model.rds"))
saveRDS(predict_conflict_xgb, file.path(OUT_XGB, "predict_conflict_fn.rds"))

cat("\n=== XGBoost 模型完成！结果保存在:", OUT_XGB, "===\n")

# ==================== Random Forest 模型部分 ====================
cat("\n=== Training Random Forest Model ===\n")

# 13. 预处理数据（为randomForest准备）
prepped <- prep(rec, training = train_df)
train_processed <- bake(prepped, train_df) %>%
  mutate(CLASS = as.factor(CLASS))
test_processed <- bake(prepped, test_df) %>%
  mutate(CLASS = as.factor(CLASS))

# 14. 训练随机森林模型（使用经验参数）
# 计算特征数量
n_features <- ncol(train_processed) - 1

# 设置经验参数
best_params <- list(
  mtry = floor(sqrt(n_features)),  # 常用规则：特征数的平方根
  ntree = 500,                     # 500棵树通常足够
  nodesize = 5                     # 终端节点最小样本数
)

cat(sprintf("Using empirical parameters: mtry=%d, ntree=%d, nodesize=%d\n",
            best_params$mtry, best_params$ntree, best_params$nodesize))

# 准备训练数据
x_train <- train_processed %>% select(-CLASS)
y_train <- train_processed$CLASS
x_test <- test_processed %>% select(-CLASS)
y_test <- test_processed$CLASS

# 训练最终随机森林模型
cat("Training final Random Forest model...\n")
rf_final <- randomForest(
  x = x_train,
  y = y_train,
  ntree = best_params$ntree,
  mtry = best_params$mtry,
  nodesize = best_params$nodesize,
  importance = TRUE,
  do.trace = 50,  # Show progress every 50 trees
  proximity = FALSE,
  keep.forest = TRUE
)

# 15. 测试集性能评估
pred_prob_rf <- predict(rf_final, x_test, type = "prob")[, "1"]

roc_obj_rf <- pROC::roc(y_test, pred_prob_rf)
auc_val_rf <- as.numeric(roc_obj_rf$auc)
pr_obj_rf  <- PRROC::pr.curve(scores.class0 = pred_prob_rf[y_test == "1"],
                              scores.class1 = pred_prob_rf[y_test == "0"], curve = TRUE)

# 保存性能指标
writeLines(
  sprintf("Best Parameters: mtry=%d, ntree=%d, nodesize=%d\nTest AUROC: %.4f\nTest AUPRC: %.4f",
          best_params$mtry, best_params$ntree, best_params$nodesize, 
          auc_val_rf, pr_obj_rf$auc.integral),
  file.path(OUT_RF, "performance.txt")
)

cat(sprintf("\nRandom Forest Test Set Performance:\nAUROC: %.4f\nAUPRC: %.4f\n", auc_val_rf, pr_obj_rf$auc.integral))

# 16. 特征重要性分析
importance_df <- as.data.frame(rf_final$importance) %>%
  mutate(feature = rownames(.)) %>%
  arrange(desc(MeanDecreaseGini))

# 取前20个重要特征
top_features <- head(importance_df, 20)

# 绘制特征重要性图
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

print(p_importance)
ggsave(file.path(OUT_RF, "feature_importance.png"), p_importance, width = 10, height = 8, dpi = 300)

# 保存特征重要性结果
write.csv(importance_df, file.path(OUT_RF, "feature_importance.csv"), row.names = FALSE)

# 17. 误差曲线（OOB误差随树数量变化）
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
  
  print(p_error)
  ggsave(file.path(OUT_RF, "error_rates.png"), p_error, width = 10, height = 6, dpi = 300)
}

# 18. 混淆矩阵
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

print(p_conf)
ggsave(file.path(OUT_RF, "confusion_matrix.png"), p_conf, width = 8, height = 6, dpi = 300)

# 计算评估指标
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
sensitivity <- conf_matrix[2, 2] / sum(conf_matrix[, 2])  # Recall/Sensitivity
specificity <- conf_matrix[1, 1] / sum(conf_matrix[, 1])  # Specificity
precision <- conf_matrix[2, 2] / sum(conf_matrix[2, ])    # Precision

metrics_df <- data.frame(
  Metric = c("Accuracy", "Sensitivity (Recall)", "Specificity", "Precision", "AUROC", "AUPRC"),
  Value = c(accuracy, sensitivity, specificity, precision, auc_val_rf, pr_obj_rf$auc.integral)
)

write.csv(metrics_df, file.path(OUT_RF, "evaluation_metrics.csv"), row.names = FALSE)

# 19. 封装预测函数
predict_conflict_rf <- function(newdata){
  # Preprocess new data
  newdata_processed <- bake(prepped, newdata)
  
  # Remove CLASS column if present
  if("CLASS" %in% colnames(newdata_processed)) {
    newdata_processed <- newdata_processed %>% select(-CLASS)
  }
  
  # Predict probabilities
  pred_probs <- predict(rf_final, newdata_processed, type = "prob")
  
  # Return probability of class 1
  if(ncol(pred_probs) == 2 && "1" %in% colnames(pred_probs)) {
    return(pred_probs[, "1"])
  } else if(ncol(pred_probs) == 2) {
    return(pred_probs[, 2])  # Assume second column is class 1
  } else {
    return(pred_probs)
  }
}

# 20. 保存模型和预处理对象
saveRDS(rf_final, file.path(OUT_RF, "final_rf_model.rds"))
saveRDS(prepped, file.path(OUT_RF, "preprocessing_recipe.rds"))
saveRDS(predict_conflict_rf, file.path(OUT_RF, "predict_conflict_rf_fn.rds"))

# 21. 输出总结
cat("\n=== Random Forest Training Complete ===\n")
cat("Models saved to:", OUT_RF, "\n")
cat("Files generated:\n")
cat("1. final_rf_model.rds - Trained Random Forest model\n")
cat("2. preprocessing_recipe.rds - Data preprocessing recipe\n")
cat("3. predict_conflict_rf_fn.rds - Prediction function\n")
cat("4. performance.txt - Performance metrics\n")
cat("5. feature_importance.csv/png - Feature importance\n")
cat("6. evaluation_metrics.csv - Detailed evaluation metrics\n")
cat("7. confusion_matrix.png - Confusion matrix\n")

if(exists("err_rate_df")) {
  cat("8. error_rates.png - Error rate curves\n")
}

cat("\n=== 全部完成！ ===\n")
cat("XGBoost 结果保存在:", OUT_XGB, "\n")
cat("Random Forest 结果保存在:", OUT_RF, "\n")
cat("\n使用示例:\n")
cat("1. XGBoost模型: model_xgb <- readRDS('output/ml/final_xgb_model.rds')\n")
cat("2. Random Forest模型: model_rf <- readRDS('output/rf/final_rf_model.rds')\n")
cat("3. 预测新数据: pred_xgb <- predict_conflict_xgb(new_data)\n")
cat("4. 预测新数据: pred_rf <- predict_conflict_rf(new_data)\n")