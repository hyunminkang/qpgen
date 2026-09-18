# qpgentools match-prs-pheno

## Summary 

`qpgentools match-prs-pheno` matches samples between PRS (Polygenic Risk Score) and phenotype matrices. It is typically used for sample identity verification or to detect sample swaps by correlating genetically predicted risks (PRS) with observed phenotypes across multiple traits. The tool identifies the best matching PRS profile for each phenotype sample.

A typical running example command is given below:

```bash
qpgentools match-prs-pheno --prs [prs_matrix] --pheno [pheno_matrix] --out [out_prefix]
```

## Required options

* `--prs` : Input PRS matrix file. The file should contain PRS values for multiple traits across samples.
* `--pheno` : Input phenotype matrix file. The file should contain observed phenotype values for the corresponding traits.
* `--out` : Output prefix to store the matching results.

## Additional Options

* `--cov` : Input covariate matrix file (optional). If provided, phenotypes will be adjusted for covariates before matching.
* `--sample-tsv` : TSV file containing sample ID mapping between PRS and phenotype files. Format: `[PRS_SAMPLE_ID] [PHENO_SAMPLE_ID]`. 
* `--trait-tsv` : TSV file containing trait ID mapping between PRS and phenotype files. Format: `[PRS_TRAIT_ID] [PHENO_TRAIT_ID]`.
* `--weights` : Input weights file for each phenotype. Format: `[PHENO_ID] [WEIGHT]`. If not provided, weights are estimated based on the correlation between PRS and phenotype in overlapping samples.
* `--prs-format` : Format of the PRS file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'.
* `--pheno-format` : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'.
* `--cov-format` : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'.
* `--missing-str` : Comma-separated strings that represent missing values in the phenotype and covariate matrices (default: `NA`). The PRS matrix must not contain missing values; the program stops with an error if it does.
* `--cov-impute-mean` : When covariates contain missing values, mean-impute them (using the mean of observed samples per covariate) instead of dropping samples with any missing covariate (default: false). See [Handling missing values](#handling-missing-values).
* `--rint` : Perform rank-based inverse normal transformation on phenotypes after covariate adjustment (default: false).
* `--mahalanobis` : Use Mahalanobis distance for matching, which accounts for correlations between traits (default: false).
* `--weight-prs-mh` : Weight given to the PRS covariance when the PRS and phenotype covariances are blended into the Mahalanobis metric, as `w * COV_prs + (1 - w) * COV_pheno` (default: 0.5). Only used with `--mahalanobis`.
* `--ledoit-wolf` : Estimate the shrinkage intensity for the covariance matrix by Ledoit-Wolf shrinkage instead of using a fixed `--lambda` (default: false). Requires `--lambda 0`; combining a non-zero `--lambda` with `--ledoit-wolf` is an error. Only used with `--mahalanobis`.
* `--lambda` : Fixed shrinkage parameter for the Mahalanobis covariance, applied as `(1 - lambda) * COV + lambda * I`. Must be between 0 and 1 (default: 0.0). Only used with `--mahalanobis`.
* `--mh-exact-norm` : With `--mahalanobis`, compute the PRS norm exactly for each missingness pattern found among the phenotyped samples, instead of approximating it with all traits (default: false). Exact norms cost one `n_prs x p^2` matrix product per unique missingness pattern, so this can be slow when many patterns exist. Only relevant when the phenotype matrix has missing values.
* `--min-weight` : Minimum weight (correlation) per trait to include in the analysis; traits below this are given zero weight (default: 0.0).
* `--z-threshold` : Z-score threshold for declaring a lenient match (default: 1.96).
* `--z-diff` : Z-score difference between the best and second-best match to declare a clear match (default: 2.0).
* `--threads` : Number of threads used by Eigen for the matrix operations (default: 1).

## Handling missing values

Missing phenotype values (strings listed in `--missing-str`, `NA` by default) are ignored rather than imputed. Concretely:

* Each trait is standardized using the mean and standard deviation of its observed values only. Missing cells contribute nothing to any downstream inner product.
* Per-trait weights (`Weight` in the weights file) are estimated from the matched samples in which the phenotype is observed; `N.Obs` reports that count.
* For each phenotyped sample, the weighted correlation is computed over its observed traits only, and normalized by the sum of absolute weights of those traits. Z-scores, ranks and match statuses therefore equal what would be obtained by dropping the missing traits for that sample. `N.Traits` reports the number of observed traits with non-zero weight.
* With `--mahalanobis`, the trait covariance is estimated from the standardized matrix with missing cells at zero (mean imputation). The cross term between a PRS sample and a phenotyped sample is restricted to the observed traits on both sides. The PRS norm is approximated using all traits by default; `--mh-exact-norm` computes it exactly per missingness pattern.
* With `--cov`, traits with missing values are adjusted by regressing on the covariates using only their observed samples. Traits without missing values are adjusted together in a single batch.
* With `--rint`, the rank-based inverse normal transformation is applied to the observed values of each trait.
* A phenotyped sample with no observed trait of non-zero weight is reported with `MatchStatus` `NO_OBS_TRAITS`, `N.Traits` 0 and `NA` in all numeric columns.

Missing covariate values are handled by dropping every sample with at least one missing covariate from the phenotype matrix before adjustment. The number of dropped samples is reported in the log, and these samples do not appear in the output. Use `--cov-impute-mean` to keep them with mean-imputed covariates instead.

The PRS matrix is assumed to be complete. If it contains a value matching `--missing-str`, the program stops with an error.

There is currently no minimum number of observed samples per trait or observed traits per sample; adding such thresholds is planned.

## Expected Output

The following files are expected to be generated:

* `[out_prefix].weights.tsv.gz` : (If `--weights` is not provided) Calculated weights for each trait.
    * Columns: `Trait`, `Weight`, `N.Obs` (number of matched samples with observed phenotype used to estimate the weight)
* `[out_prefix].match.assigned.tsv.gz` : Summary of the best match for each phenotype sample.
    * Columns:
        * `ID.Pheno`: Sample ID in the phenotype file.
        * `N.Traits`: Number of observed traits with non-zero weight used for this sample.
        * `MatchStatus`: Status of the match (e.g., `BEST_MATCH`, `NO_MATCH`, `LENIENT_MATCH`, or `NO_OBS_TRAITS` when no trait is observed).
        * `ID.self`: Corresponding sample ID in the PRS file (based on mapping or ID match).
        * `Z.self`: Z-score of the match with itself.
        * `COR.self`: Weighted correlation with itself.
        * `Rank.self`: Rank of the self-match among all PRS samples.
* `[out_prefix].match.all.tsv.gz` : Detailed output including the top 5 matches for each phenotype sample.
    * Columns: `ID.Pheno`, `N.Traits`, `MatchStatus`, `ID.self`, `Z.self`, `COR.self`, `Rank.self`, followed by `ID.1st`, `Z.1st`, `COR.1st` through `ID.5th`, `Z.5th`, `COR.5th` for the five best-matching PRS samples.

## Full Usage 

The full usage of `qpgentools match-prs-pheno` can be viewed with the `--help` option:

```
$ ./qpgentools match-prs-pheno --help
[bin/qpgentools match-prs-pheno] -- Match PRS and phenotype matrices based on overlapping samples and compute weights for each phenotype

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input Options ==
   --prs             [STR: ]             : Input PRS file
   --pheno           [STR: ]             : Input phenotype matrix
   --cov             [STR: ]             : Input covariate matrix (optional)
   --sample-tsv      [STR: ]             : TSV file containing input sample IDs PRS and phenotype files in [PRS_SAMPLE_ID] [PHENO_SAMPLE_ID] format. If they use the same IDs, use only a single column if subsetting samples are needed
   --trait-tsv       [STR: ]             : TSV file containing input trait IDs PRS and phenotype files in [PRS_TRAIT_ID] [PHENO_TRAIT_ID] format. If they use the same IDs, use only a single column if subsetting traits are needed
   --weights         [STR: ]             : Input weights file for each phenotype in [PHENO_ID] [WEIGHT] format
   --prs-format      [STR: regenie]      : Format of the PRS file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'
   --pheno-format    [STR: regenie]      : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
   --cov-format      [STR: regenie]      : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'
   --missing-str     [STR: NA]           : Comma-separated strings representing missing values in the phenotype and covariate matrices (default: 'NA'). The PRS matrix must not contain missing values

== Output options ==
   --out             [STR: ]             : Output prefix

== Analysis options ==
   --cov-impute-mean [FLG: OFF]          : Mean-impute missing covariate values instead of dropping samples with any missing covariate (default: false)
   --rint            [FLG: OFF]          : Perform rank-based inverse normal transformation after covariate adjustment (default: false)
   --mahalanobis     [FLG: OFF]          : Use Mahalanobis distance for matching (default: false)
   --ledoit-wolf     [FLG: OFF]          : Use Ledoit-Wolf shrinkage for covariance estimation when using Mahalanobis distance (default: false)
   --mh-exact-norm   [FLG: OFF]          : With --mahalanobis, compute PRS norms exactly for each phenotype missingness pattern instead of approximating with all traits (slower when many patterns exist; default: false)
   --weight-prs-mh   [FLT: 0.50]         : Weight for PRS distance when combining with weighted correlation (default: 0.5)
   --lambda          [FLT: 0.00]         : Regularization parameter for Mahalanobis distance, between 0 and 1 (default: 0.0)
   --min-weight      [FLT: 0.00]         : Minimum weight (in r) per trait to set to zero (default: 0.0)
   --z-threshold     [FLT: 1.96]         : Z-score threshold for lenient matching (default: 1.96)
   --z-diff          [FLT: 2.00]         : Z-score difference to declare a clear match (default: 2.0)
   --threads         [INT: 1]            : Number of threads to use (default: 1)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```