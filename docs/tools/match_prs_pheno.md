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
* `--rint` : Perform rank-based inverse normal transformation on phenotypes after covariate adjustment (default: false).
* `--mahalanobis` : Use Mahalanobis distance for matching, which accounts for correlations between traits (default: false).
* `--weight-prs-mh` : Weight given to the PRS covariance when the PRS and phenotype covariances are blended into the Mahalanobis metric, as `w * COV_prs + (1 - w) * COV_pheno` (default: 0.5). Only used with `--mahalanobis`.
* `--ledoit-wolf` : Estimate the shrinkage intensity for the covariance matrix by Ledoit-Wolf shrinkage instead of using a fixed `--lambda` (default: false). Requires `--lambda 0`; combining a non-zero `--lambda` with `--ledoit-wolf` is an error. Only used with `--mahalanobis`.
* `--lambda` : Fixed shrinkage parameter for the Mahalanobis covariance, applied as `(1 - lambda) * COV + lambda * I`. Must be between 0 and 1 (default: 0.0). Only used with `--mahalanobis`.
* `--min-weight` : Minimum weight (correlation) per trait to include in the analysis; traits below this are given zero weight (default: 0.0).
* `--z-threshold` : Z-score threshold for declaring a lenient match (default: 1.96).
* `--z-diff` : Z-score difference between the best and second-best match to declare a clear match (default: 2.0).
* `--threads` : Number of threads used by Eigen for the matrix operations (default: 1).

## Expected Output

The following files are expected to be generated:

* `[out_prefix].weights.tsv.gz` : (If `--weights` is not provided) Calculated weights for each trait.
    * Columns: `Trait`, `Weight`
* `[out_prefix].match.assigned.tsv.gz` : Summary of the best match for each phenotype sample.
    * Columns:
        * `ID.Pheno`: Sample ID in the phenotype file.
        * `MatchStatus`: Status of the match (e.g., `BEST_MATCH`, `NO_MATCH`, `LENIENT_MATCH`).
        * `ID.self`: Corresponding sample ID in the PRS file (based on mapping or ID match).
        * `Z.self`: Z-score of the match with itself.
        * `COR.self`: Weighted correlation with itself.
        * `Rank.self`: Rank of the self-match among all PRS samples.
* `[out_prefix].match.all.tsv.gz` : Detailed output including the top 5 matches for each phenotype sample.
    * Columns: `ID.Pheno`, `MatchStatus`, `ID.self`, `Z.self`, `COR.self`, `Rank.self`, followed by `ID.1st`, `Z.1st`, `COR.1st` through `ID.5th`, `Z.5th`, `COR.5th` for the five best-matching PRS samples.

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
   --prs           [STR: ]             : Input PRS file
   --pheno         [STR: ]             : Input phenotype matrix
   --cov           [STR: ]             : Input covariate matrix (optional)
   --sample-tsv    [STR: ]             : TSV file containing input sample IDs PRS and phenotype files in [PRS_SAMPLE_ID] [PHENO_SAMPLE_ID] format. If they use the same IDs, use only a single column if subsetting samples are needed
   --trait-tsv     [STR: ]             : TSV file containing input trait IDs PRS and phenotype files in [PRS_TRAIT_ID] [PHENO_TRAIT_ID] format. If they use the same IDs, use only a single column if subsetting traits are needed
   --weights       [STR: ]             : Input weights file for each phenotype in [PHENO_ID] [WEIGHT] format
   --prs-format    [STR: regenie]      : Format of the PRS file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'
   --pheno-format  [STR: regenie]      : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
   --cov-format    [STR: regenie]      : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'

== Output options ==
   --out           [STR: ]             : Output prefix

== Analysis options ==
   --rint          [FLG: OFF]          : Perform rank-based inverse normal transformation after covariate adjustment (default: false)
   --mahalanobis   [FLG: OFF]          : Use Mahalanobis distance for matching (default: false)
   --ledoit-wolf   [FLG: OFF]          : Use Ledoit-Wolf shrinkage for covariance estimation when using Mahalanobis distance (default: false)
   --weight-prs-mh [FLT: 0.50]         : Weight for PRS distance when combining with weighted correlation (default: 0.5)
   --lambda        [FLT: 0.00]         : Regularization parameter for Mahalanobis distance, between 0 and 1 (default: 0.0)
   --min-weight    [FLT: 0.00]         : Minimum weight (in r) per trait to set to zero (default: 0.0)
   --z-threshold   [FLT: 1.96]         : Z-score threshold for lenient matching (default: 1.96)
   --z-diff        [FLT: 2.00]         : Z-score difference to declare a clear match (default: 2.0)
   --threads       [INT: 1]            : Number of threads to use (default: 1)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```