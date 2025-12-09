#' Performance Metrics of Candidate Models
#'
#' A dataset containing key performance indicators (e.g., AUROC, Sensitivity)
#' for the candidate models evaluated during package development.
#'
#' @format A data frame with 7 rows and 4 columns:
#' \describe{
#'   \item{Indicator}{The name of the performance metric.}
#'   \item{XGBoost}{Performance values for the XGBoost model.}
#'   \item{RF (90% Training)}{Performance values for the Random Forest model trained on 90% of the data.}
#'   \item{RF (80% Training)}{Performance values for the Random Forest model trained on 80% of the data.}
#' }
"performance_metrics"

#' Sample Variant Data for Pathogenicity Prediction
#'
#' A subset of the ClinVar Conflicting dataset containing 50 sample variants.
#' This data is used for examples and demonstrations in the package.
#'
#' @format A data frame with 50 rows and 45 columns, including key features
#'   like CADD_PHRED, CADD_RAW, and the true 'CLASS' label.
"sample_variants"
