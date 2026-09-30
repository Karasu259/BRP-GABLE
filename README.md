# BRP-GABLE

This repository contains the R code, numerical results, and figures accompanying the study of Bootstrap Rule Pooling for Generalized ABLE (BRP-GABLE). BRP-GABLE pools rules generated from bootstrap-specific GABLE fits and performs a single global lasso selection over the pooled rules, with or without truncated linear terms.

The code covers binary, continuous, and time-to-event outcomes. It includes the main simulation study, supplementary simulations, and both MGUS2 real-data experiments.

## Preprint

The accompanying manuscript is available on arXiv:

Li Y, Wan K, Shimokawa T, and Tanioka K.  
BRP-GABLE: An Interpretable Rule-Based Prediction Framework for Multiple Outcome Types.  
arXiv:2609.36839, 2026.  
https://arxiv.org/abs/2609.36839

## Repository structure

```text
R/                         Core rule-generation, model-fitting, and analysis functions
scripts/simulation/        Main simulations with n = 1,000
scripts/appendix/          Sample-size and RuleFit-depth supplementary simulations
scripts/real_data/         MGUS2 Experiments 1 and 2
scripts/plots/             Scripts used to reproduce all released figures
results/simulation/        Main and supplementary simulation results
results/real_data/         MGUS2 results
figures/simulation/        Simulation figures used in the manuscript and appendix
figures/real_data/         MGUS2 figures used in the manuscript and appendix
```

The run scripts automatically locate and switch to the repository root when executed with `Rscript`, sourced in RStudio, or run line by line in RStudio. The plotting commands below should be run from the repository root.

## Software requirements

The analyses were run with R 4.4.2 on Windows 11. Required R packages are:

```r
install.packages(c(
  "survival", "glmnet", "MASS", "pre", "pROC", "gbm",
  "randomForest", "randomForestSRC", "foreach", "doParallel",
  "doRNG", "ggplot2", "patchwork", "scales", "rstudioapi"
))
```

The baseline implementations used `pre` 1.0.7, `gbm` 2.3.1, `randomForest` 4.7-1.2, and `randomForestSRC` 3.6.2. The result manifests under `results/simulation/` contain the recorded R session information for the supplementary simulation runs.

## Main simulation study

The main study used four scenarios for each outcome type, 100 replications per scenario, a total sample size of 1,000, and a 1:1 training--test split. The fixed BRP-GABLE settings were 50 bootstrap replications, a maximum of 10 generated rule terms per bootstrap fit, maximum rule degree 2, minimum child-rule support of 5%, and a 20%-of-sample-size spacing for the admissible rank grid. The final lasso penalty was selected by 10-fold cross-validation using the one-standard-error rule.

The seven evaluated models were BRP-GABLE full, BRP-GABLE rule-only, GABLE, RuleFit, the outcome-specific linear model, gradient boosting machines, and random forests.

Run one outcome at a time:

```sh
Rscript scripts/simulation/run_binary.R
Rscript scripts/simulation/run_continuous.R
Rscript scripts/simulation/run_survival.R
```

Each script is configured for 20 parallel workers. Adjust `n_workers` in the corresponding script when fewer cores are available. Full reproduction is computationally intensive.

Recreate the three main figures from the included CSV files:

```sh
Rscript scripts/plots/plot_formal_simulation_results.R
```

## Supplementary simulations

The sample-size analyses change only the total sample size from the main simulation design.

```sh
Rscript scripts/appendix/run_appendix_n500.R --workers=20
Rscript scripts/appendix/run_appendix_n2000.R --workers=20
```

A short validation run is available by adding `--smoke` to either command.

The matched-depth analysis refits only RuleFit with maximum tree depth 2 while retaining the data, splits, and seeds used in the main simulations:

```sh
Rscript scripts/appendix/13_run_rulefit_maxdepth2_appendix.R --workers=20
```

The supplementary figures can be recreated with:

```sh
Rscript scripts/plots/plot_appendix_n500_results.R
Rscript scripts/plots/plot_appendix_n2000_results.R
Rscript scripts/plots/plot_appendix_sample_size_comparison.R
Rscript scripts/plots/plot_appendix_rulefit_maxdepth2.R
```

## MGUS2 real-data analysis

The data are obtained directly from `survival::mgus2`; no private data are required. The analysis uses overall survival (`futime` and `death`) and the predictors age, sex, year of diagnosis, hemoglobin, serum creatinine, and M-spike. Complete-case preprocessing gives 1,338 observations.

Experiment 1 performs 100 repeated random 7:3 training--test splits and evaluates the same seven models as the simulation study:

```sh
Rscript scripts/real_data/run_real_survival.R
Rscript scripts/plots/plot_real_data_experiment1.R
```

Experiment 2 fits BRP-GABLE once to the complete dataset and produces the retained-term, rule support--effect, patient activation, rule similarity, ALE, and Kaplan--Meier outputs:

```sh
Rscript scripts/real_data/run_real_survival_full.R
Rscript scripts/plots/plot_real_data_experiment2.R
```

The supplementary ALE comparison with RuleFit, GBM, and random forests is generated with:

```sh
Rscript scripts/plots/plot_appendix_ale_rulefit_comparison.R
```

This last script uses the included BRP-GABLE fit and refits missing baseline-model caches when necessary.

## Included results

The repository includes the CSV outputs needed for the released figures, together with PDF and PNG versions of the manuscript figures. Large simulation-model objects and cached baseline fits are intentionally excluded. The small full-data BRP-GABLE fit required by the Experiment 2 plotting code is retained under `results/real_data/experiment2/`.

The principal result columns are:

- `metric_value`: AUC for binary outcomes, mean squared error for continuous outcomes, or C-index for time-to-event outcomes.
- `n_rules_total`: total number of selected final terms, including retained linear terms for rule-based methods with an explicit linear component.
- `n_rule_terms_selected` and `n_linear_selected`: selected rule and linear terms separately.
- `prop_deg1`: proportion of degree-1 terms among the selected final terms.
- `avg_rule_degree`: average degree of the selected final terms.

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.
