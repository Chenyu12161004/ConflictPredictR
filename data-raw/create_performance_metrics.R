## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

indicators <- c("AUROC", "AUPRC", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1-Score")

xgboost_metrics <- c(0.7382, 0.4571, 0.7520, 0.2716, 0.9139, 0.5153, 0.2200)
rf_90_metrics <- c(0.7397, 0.4612, 0.7530, 0.2944, 0.9076, 0.5179, 0.2355)
rf_80_metrics <- c(0.7514, 0.4575, 0.6999, 0.6501, 0.7166, 0.4361, 0.5233)

performance_metrics <- data.frame(
  Indicator = indicators,
  XGBoost = xgboost_metrics,
  `RF (90% Training)` = rf_90_metrics,
  `RF (80% Training)` = rf_80_metrics,
  check.names = FALSE
)

usethis::use_data(performance_metrics, overwrite = TRUE)
