README
================

# PathoPrediction: A Focused Pathogenicity Predictor

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/Chenyu12161004/PathoPrediction/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Chenyu12161004/PathoPrediction?branch=main)
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Overview

`PathoPrediction` is a tool for predicting the pathogenicity of human
genetic variants based on a set of genomic and computational features,
using pre-trained Random Forest model(at 80% split). This model was
chosen over other alternatives, XGBoost and Random Forest (at 90%
split), for its higher AUROC and especially better sensitivity in
identifying potentially harmful variants. The package provides functions
for prediction. It also offers visualizations to explain model behavior,
including feature importance bar plots and diagnostic charts (OOB error
and confusion matrix).

## Key Features

- **Focused Prediction**: Deliver pathogenicity predictions using a
  single, robust Random Forest model. This approach prioritizes high
  sensitivity to ensure potentially harmful variants are effectively
  identified.

- **Clear Feature Insights**: Help to understand model behavior with
  clear plots for feature importance (`plot_rf_importance()`) and model
  diagnostics (`plot_rf_diagnostics()`).

- **User-Friendly Data Validation**: Include a built-in validation
  function, `validate_input_data()`, that checks the input data for
  common formatting issues and provides clear error messages.

- **Rigorously Tested**: The package is equipped with a suite of
  automated tests to ensure functions behave reliably and predictably,
  giving users confidence in the results.

- **Transparent Model Selection**: The package includes a performance
  comparison table, openly showing why the Random Forest model (at 80%
  split) was chosen over XGBoost and Random Forest (at 90% split) for
  this specific task.

## Model Selection

This package is built around a single Random Forest (RF) model (at 80%
split). During development, other models like XGBoost were also
evaluated. While all these models showed well enough overall
performance, the Random Forest model demonstrated higher Sensitivity.
This means it is better at correctly identifying true pathogenic
variants, which minimizes the risk of missing a potentially significant
variant, providing a more cautious and reliable result.

    ## ℹ Loading PathoPrediction

    ## Warning: package 'testthat' was built under R version 4.5.2

| Indicator   | XGBoost | RF (90% Training) | RF (80% Training) |
|:------------|--------:|------------------:|------------------:|
| AUROC       |  0.7382 |            0.7397 |            0.7514 |
| AUPRC       |  0.4571 |            0.4612 |            0.4575 |
| Accuracy    |  0.7520 |            0.7530 |            0.6999 |
| Sensitivity |  0.2716 |            0.2944 |            0.6501 |
| Specificity |  0.9139 |            0.9076 |            0.7166 |
| Precision   |  0.5153 |            0.5179 |            0.4361 |
| F1-Score    |  0.2200 |            0.2355 |            0.5233 |

Model Performance Comparison

------------------------------------------------------------------------

## Model Performance

The performance of the Random Forest model (at 80% split) was evaluated
using AUROC and AUPRC.

|           ROC and PRC Curves            |
|:---------------------------------------:|
| ![](man/figures/roc_and_prc_curves.png) |
|     *AUROC: 0.7514, AUPRC: 0.4575*      |

------------------------------------------------------------------------

## Installation

You can install the development version of `PathoPrediction` from GitHub
with:

``` r
# install.packages("devtools")
# devtools::install_github("Chenyu12161004/PathoPrediction")
```

## Usage

First, load the package:

``` r
devtools::load_all()
```

    ## ℹ Loading PathoPrediction

### Data validation

Before making predictions, it is a good practice to validate the input
data. The `validate_input_data()` function checks for required columns
and data types.

``` r
# Validate the built-in sample data (this should pass)
validate_input_data(sample_variants)
```

    ## Starting Data Validation

    ## All required columns are present.

    ## Warning in recipes::bake(recipe_obj, new_data = newdata): ! There were 31 columns that were factors when the recipe was prepped:
    ## • `CHROM`, `REF`, `ALT`, `CLNDISDB`, `CLNDISDBINCL`, `CLNDN`, `CLNDNINCL`,
    ##   `CLNHGVS`, `CLNSIGINCL`, `CLNVC`, `CLNVI`, `MC`, `Allele`, `Consequence`,
    ##   `IMPACT`, `SYMBOL`, `Feature_type`, `Feature`, …, `MOTIF_NAME`, and
    ##   `HIGH_INF_POS`
    ## ℹ This may cause errors when processing new data.

    ## Data is compatible with the model's preprocessing steps.

    ## Validation Successful

### Predict Pathogenicity

Use `predict_pathogenicity()` on the sample data to get predictions. The
function returns the predicted class and the probability of the variant
being pathogenic (class “1”).

``` r
# Use the sample data included in the package
predictions <- predict_pathogenicity(sample_variants)
```

    ## Loading Random Forest model and preprocessing recipe...

    ## Processing new data...

    ## Warning in recipes::bake(recipe_obj, new_data = newdata): ! There were 31 columns that were factors when the recipe was prepped:
    ## • `CHROM`, `REF`, `ALT`, `CLNDISDB`, `CLNDISDBINCL`, `CLNDN`, `CLNDNINCL`,
    ##   `CLNHGVS`, `CLNSIGINCL`, `CLNVC`, `CLNVI`, `MC`, `Allele`, `Consequence`,
    ##   `IMPACT`, `SYMBOL`, `Feature_type`, `Feature`, …, `MOTIF_NAME`, and
    ##   `HIGH_INF_POS`
    ## ℹ This may cause errors when processing new data.

    ## Making predictions...

    ## Done.

``` r
# View the key prediction results
result_subset <- predictions[, c("CADD_PHRED", "CADD_RAW", ".pred_class", ".pred_prob_pathogenic")]
knitr::kable(head(result_subset))
```

| CADD_PHRED | CADD_RAW | .pred_class | .pred_prob_pathogenic |
|-----------:|---------:|:------------|----------------------:|
|     17.630 | 2.218361 | 1           |                 0.666 |
|     23.300 | 3.717784 | 0           |                 0.014 |
|     16.280 | 2.010398 | 1           |                 0.520 |
|     28.000 | 6.029584 | 0           |                 0.196 |
|     25.300 | 5.093486 | 0           |                 0.374 |
|      8.856 | 0.704474 | 1           |                 0.714 |

Explanation of the new columns in the report:

- `.pred_class`: The final classification from the Random Forest model
  (1 = pathogenic, 0 = benign).
- `pred_prob_pathogenic`: The probability that the variant is pathogenic
  (class “1”), according to the model.

### Visualization: Understand Model Behavior

The package provides two functions to visualize what drives the
decisions of the model.

#### Feature Importance

`plot_rf_importance()` generates a bar chart showing which features had
the most influence on the model’s predictions.

``` r
plot_rf_importance()
```

    ## Loading Random Forest model...

![](README_files/figure-gfm/importance-plot-1.png)<!-- -->

#### Model Diagnostics

`plot_rf_diagnostics()` demonstrates the internal performance of the
model, showing the Out-of-Bag (OOB) error rate during training and the
final confusion matrix.

- **Out-of-Bag (OOB) Error Rate**: This plot shows the error rate of the
  model as it was being trained. It is generated from data stored within
  the model object itself and is therefore static. It reflects the
  internal performance of the model during its creation.

- **Confusion Matrix**: This plot evaluates performance of the model on
  a specific dataset. By default, it uses the built-in sample_variants
  data. If you provide your own newdata, this plot will update to show
  the accuracy of the model on that specific data.

``` r
plot_rf_diagnostics()
```

    ## Loading Random Forest model...

    ## No 'newdata' provided. Displaying confusion matrix from OOB samples.

![](README_files/figure-gfm/diagnostics-plot-1.png)<!-- -->

------------------------------------------------------------------------

## About the Model

This package is built around a single machine learning model trained on
the ClinVar Conflicting dataset.

- **Data Preprocessing**: All data was prepared using a `recipes`
  pipeline. This process included handling missing values, converting
  categorical features into dummy variables, removing zero-variance
  predictors, and reducing high correlation between numeric features.

- **Random Forest Model**: Trained using the `randomForest` package with
  standard, empirically chosen parameters (`mtry` = sqrt(p), `ntree` =
  500).

The functions in this package are direct adaptations of the original
analysis script, ensuring that the research workflow is clear and
reproducible.
