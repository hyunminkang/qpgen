# qpgentools pair-assoc

## Summary 

`qpgentools pair-assoc` performs rapid association tests between selected variant-trait pairs.

A typical running example command is given below:

```bash
qpgentools pair-assoc --pgen-list [list] --pheno [pheno] --cov [cov] --pairs [pairs] --out [out_prefix]
```

## Required options

* Input genotype files : Either `--pgen-list` or `--pgen`, `--psam`, `--pivar` options are required. See [Genotype file formats](../formats/genotypes.md) for more details.
* `--pheno` : Input phenotype matrix in Regenie, TensorQTL, or TSV format. See [Phenotype file formats](../formats/phenotypes.md) for more details.
* `--pairs` : Input TSV file containing `[trait_id]` and `[variant_id]` or `[trait_id]` and `[region]` at each line to perform association tests
* `--out` : Output TSV file to store the association analysis output 

## Additional Options

* `--cov` : Input covariate matrix in Regenie, TensorQTL, or TSV format. See [Phenotype file formats](../formats/phenotypes.md) for more details. When provided, **both** the phenotypes and the genotypes are residualized against the covariates (see [Covariate handling](#covariate-handling)).
* `--sample` : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
* `--pheno-format` : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
* `--cov-format` : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
* `--skip-rint` : Skip tests based on rank-based inverse normal transformation (default: false). With `--cov`, the RINT is applied to the covariate-adjusted phenotype, as `rect-assoc --rint-after-adj` does.
* `--offset-pheno` : The number of index columns (i.e. not containing values) in the phenotype files (default: 1)
* `--offset-cov` : The number of index columns (i.e. not containing values) in the covariate files (default: 1)
* `--icol-pivar-idx` : 1-based column index for the variant index in the pvar file (default: 9)
* `--icol-pheno-id` : 1-based column index for the phenotype ID in the phenotype file (default: 1)
* `--icol-cov-id` : 1-based column index for the covariate ID in the covariate file (default: 1)
* `--colname-pheno-sample` : When `--sample` is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')
* `--colname-geno-sample` : When `--sample` is provided, the column name for the sample IDs in the genotype file (default: 'geno')
* `--colname-pair-trait` : Column name of the trait ID in the pair file (default: 'trait')
* `--colname-pair-variant` : Column name of the variant ID in the pair file (default: 'variant')
* `--colname-pair-region` : Column name of the region string in the pair file (default: 'region')

## Covariate handling

When `--cov` is supplied, the phenotypes are adjusted for the covariates by linear regression
and the genotypes are residualized against the *same* covariates (Frisch-Waugh-Lovell), so
that `BETA` is the partial effect of the variant. Adjusting only the phenotype -- as earlier
versions of this command did -- leaves genotype variance that is collinear with the covariates
(e.g. genotype PCs) in the design, which shrinks `BETA` and the test statistic.

* The `_RINT` statistics use the rank-based inverse normal transform of the covariate-adjusted
  phenotype, tested against the same residualized genotypes. They match
  [`rect-assoc`](rect_assoc.md) with `--rint-after-adj`. Because the transformed phenotype is no
  longer exactly orthogonal to the covariates, its residual sum of squares on the intercept and
  covariates replaces `y'y` in the residual variance, so the `_RINT` statistics are those of the
  OLS fit `RINT(adjusted y) ~ 1 + covariates + genotype`.
* Missing hard-call genotypes are mean-imputed before residualization, as in `rect-assoc` and
  [`region-assoc`](region_assoc.md), so every sample enters the regression. `N`, `A1FREQ`,
  `N_RR`/`N_RA`/`N_AA`, and `INFO` are computed from the raw genotypes, so `N` is the number of
  samples with a called genotype.
* `SE`, `TSTAT` and `LOG10P` (and the `_RINT` columns) use the residual degrees of freedom of the
  full model `y ~ 1 + covariates + genotype`: `df = n - 2 - k`, where `n` is the number of analyzed
  samples and `k` is the number of linearly independent covariates (the numerical rank of the
  intercept and covariates, minus 1; `k` equals the number of covariates when they are full
  rank). `BETA`, `SE`, `TSTAT` and `LOG10P` therefore equal those of an ordinary least-squares
  fit with the covariates in the model. Because missing genotypes are imputed, `n` can exceed the
  `N` column for variants with missing genotypes. Without `--cov`, `df = n - 2` as before. The
  value of `k` is printed in the log.
* A variant whose genotype is (numerically) a linear combination of the covariates -- residual
  variance below `1e-8` of its total variance, e.g. when the variant itself is a covariate -- cannot
  be tested; its `BETA`, `SE`, `TSTAT`, and `LOG10P` (and `_RINT`) are reported as NaN, and the
  number of such tests is reported at the end of the run.
* A covariate that is constant in the analyzed samples (e.g. sex when `--sample` keeps one sex) or
  a linear combination of other covariates is redundant. The adjustment uses the numerical rank of
  the intercept and covariates, so a redundant covariate is adjusted for once and costs no degree
  of freedom; the results equal those obtained without it. The command prints a warning when it
  detects redundant covariates; dropping them is still recommended.

!!! note
    Effect sizes and p-values from `--cov` runs will therefore differ from those produced by
    versions of `qpgentools` released before this fix. For the same covariates, `pair-assoc`,
    [`rect-assoc`](rect_assoc.md), and [`region-assoc`](region_assoc.md) now agree on shared
    (trait, variant) pairs. Without `--cov`, the output is unchanged.

    Earlier versions also used `df = n - 2` with `--cov`, ignoring the covariates. Their `SE` was
    too small by a factor of `sqrt((n - 2 - k) / (n - 2))` and their `LOG10P` too large (about 5% in
    `TSTAT` for 66 covariates and 729 samples). `BETA` is unchanged by the degrees-of-freedom
    correction.

!!! warning
    `pair-assoc` does not recognize missing-value tokens in the phenotype or covariate files:
    a non-numeric cell such as `NA` is read as `0`, with no error. Provide phenotype and
    covariate matrices **without missing values**.

## Expected Output

The output TSV file contains the following columns, largely resembling the Regenie output format:

* `TRAIT` : Trait name
* `CHROM` : Chromosome of the variant tested
* `GENPOS` : Base position (1-based) of the variant tested
* `ID` : Variant ID in `[CHROM]:[GENPOS]:[REF]:[ALT]` format
* `ALLELE0` : Reference allele
* `ALLELE1` : Alternative allele
* `A1FREQ`: Allele frequency of alternate allele
* `N`: Sample size. With `--cov`, the number of samples with a called genotype (see [Covariate handling](#covariate-handling)); without it, the number of samples with both a phenotype value and a called genotype
* `N_RR`: Number of RR genotypes
* `N_RA`: Number of RA genotypes
* `N_AA`: Number of AA genotypes
* `INFO`: Variant quality score (typically imputation quality)
* `BETA`: Effect size for raw phenotype association (partial effect adjusted for the covariates when `--cov` is given)
* `SE`: Standard error for raw phenotype association
* `TSTAT`: Test statistic for raw phenotype association
* `LOG10P`: Log10 p-value for raw phenotype association
* `BETA_RINT`: Effect size for rank-based inverse normal transformed phenotype
* `SE_RINT`: Standard error for rank-based inverse normal transformed phenotype
* `TSTAT_RINT`: Test statistic for rank-based inverse normal transformed phenotype
* `LOG10P_RINT`: Log10 p-value for rank-based inverse normal transformed phenotype


## Full Usage 

The full usage of `qpgentools pair-assoc` can be viewed with the `--help` option:

```
$ ./qpgentools pair-assoc --help
[bin/qpgentools pair-assoc] -- Perform pairwise association analysis for specific pairs of phenotype variant pairs

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input Options ==
   --pgen                 [STR: ]             : Input PLINK2 genotype file
   --psam                 [STR: ]             : Input PLINK2 sample file
   --pivar                [STR: ]             : Input PLINK2 index pvar file (bgzipped and tabix)
   --pgen-list            [STR: ]             : Input file containing CHROM BEG END PGEN PSAM PIVAR
   --pheno                [STR: ]             : Input phenotype matrix
   --pairs                [STR: ]             : Input file containing [trait_id] [variant_id] pairs to be tested
   --sample               [STR: ]             : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
   --cov                  [STR: ]             : Input covariate matrix (optional)
   --pheno-format         [STR: regenie]      : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
   --cov-format           [STR: regenie]      : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'

== Output options ==
   --out                  [STR: ]             : Output prefix

== Auxiliary options ==
   --jump-thres-bp        [INT: 1000000]      : Jump threshold in base pairs for the variant index (default: 1000000)
   --max-chunk-vars       [INT: 1000]         : Maximum number of variants to store at once in memory (default: 1000)
   --offset-pheno         [INT: 1]            : The number of index columns (i.e. not containing values) in the phenotype files (default: 1)
   --offset-cov           [INT: 1]            : The number of index columns (i.e. not containing values) in the covariate columns (default: 1)
   --icol-pivar-idx       [INT: 9]            : 1-based column index for the variant ID in the pvar file (default: 9)
   --icol-pheno-id        [INT: 1]            : 1-based column index for the phenotype ID in the phenotype file (default: 1)
   --icol-cov-id          [INT: 1]            : 1-based column index for the covariate ID in the phenotype file (default: 1)
   --colname-pheno-sample [STR: pheno]        : When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')
   --colname-geno-sample  [STR: geno]         : When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'geno')
   --colname-pair-trait   [STR: trait]        : Column name of the trait ID in the pait file
   --colname-pair-variant [STR: variant]      : Column name of the variant ID in the pair file
   --colname-pair-region  [STR: region]       : Column name of the region string in the pair file
   --skip-rint            [FLG: OFF]          : Skip tests based on rank-based inverse normal transformation (default: false)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```