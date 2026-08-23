# qpgentools rect-assoc

## Summary 

`qpgentools rect-assoc` performs rectangular association analysis, testing a specific list of variants against a specific list of phenotypes (or all phenotypes). This is efficient for looking up associations for a subset of variants across multiple traits, producing a "rectangular" slice of the full association summary statistics.

To test *every* variant in a genomic region instead of an explicit variant list, and to fine-map
that region with SuSiE, use [`region-assoc`](region_assoc.md).

A typical running example command is given below:

```bash
qpgentools rect-assoc --pgen [pgen] --pivar [pivar] --psam [psam] --pheno [pheno] --var-list [variants] --out [out_prefix]
```

## Required options

* Input genotype files : Either `--pgen-list` or `--pgen`, `--psam`, `--pivar` options are required. See [Genotype file formats](../formats/genotypes.md) for more details.
* `--pheno` : Input phenotype matrix in Regenie, TensorQTL, or TSV format. See [Phenotype file formats](../formats/phenotypes.md) for more details.
* `--var-list` : Input file containing the list of variants to be tested. The file can contain a single column of variant IDs (e.g. `chr:pos:ref:alt`) or 4 columns (CHROM, POS, REF, ALT).
* `--out` : Output **file** to store the association analysis output. Despite the help text calling it a prefix, the value is used as the file name as-is; use a `.gz` suffix for bgzipped output.

## Additional Options

* `--pheno-list` : Input file containing the list of phenotypes to be tested (single column with phenotype IDs). If not provided, all phenotypes in the phenotype matrix will be tested. Phenotype IDs that are absent from the phenotype matrix are silently dropped.
* `--cov` : Input covariate matrix in Regenie or TSV format. See [Covariate handling](#covariate-handling) below.
* `--sample` : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files.
* `--pheno-format` : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'.
* `--cov-format` : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'.
* `--rint-before-adj` : Perform rank-based inverse normal transformation before covariate adjustment (default: false).
* `--rint-after-adj` : Perform rank-based inverse normal transformation after covariate adjustment (default: false).
* `--min-maf` : Minimum minor allele frequency to include a variant in the analysis (default: 1e-10, i.e. only monomorphic variants are excluded).
* `--min-mac` : Minimum minor allele count to include a variant in the analysis (default: 1.0).
* `--max-chunk-vars` : Maximum number of variants held in memory at once (default: 100). Variants from `--var-list` are processed in chunks of this size.
* `--jump-thres-bp` : Jump threshold in base pairs for the variant index (default: 1000000).
* `--icol-pivar-idx` : 1-based column index for the variant index in the indexed pvar file (default: 9).
* `--icol-pheno-id` / `--icol-cov-id` : 1-based column index for the phenotype/covariate ID (default: 1).

## Covariate handling

When `--cov` is supplied, the phenotype matrix is adjusted for the covariates by linear regression
and the genotype matrix is residualized against the *same* covariates (Frisch-Waugh-Lovell), so
that `BETA` is the partial effect of the variant. Adjusting only the phenotype -- as earlier
versions of this command did -- leaves genotype variance that is collinear with the covariates
(e.g. genotype PCs) in the design, which shrinks `BETA` and the test statistic.

!!! note
    Effect sizes and p-values from `--cov` runs will therefore differ from those produced by
    versions of `qpgentools` released before this fix. [`region-assoc`](region_assoc.md) applies
    the same treatment, so the two commands now agree on shared variants.

    `LOG10P` uses `df = N - 2` and does not subtract the number of covariates.

!!! warning
    Covariate adjustment and the rank-based inverse normal transformations are currently supported
    only for phenotype and covariate matrices **without missing values**; the run aborts otherwise.

## Expected Output

The output file named by `--out` (bgzipped when the name ends in `.gz`) contains the following columns. It is in a wide format where statistics for each phenotype are appended horizontally.

* `#CHROM` : Chromosome
* `GENPOS` : Genomic Position
* `ID` : Variant ID
* `ALLELE0` : Reference Allele
* `ALLELE1` : Alternative Allele
* `A1FREQ` : Frequency of Allele 1 (AC / AN)
* `N` : Number of samples with a called genotype (`N_RR + N_RA + N_AA`)
* `N_RR` : Count of Reference/Reference genotypes
* `N_RA` : Count of Reference/Alternative genotypes
* `N_AA` : Count of Alternative/Alternative genotypes
* `TEST` : Test type (e.g., "Linear")

For each phenotype `[PHENO]`, the following columns are added:

* `BETA.[PHENO]` : Effect size
* `SE.[PHENO]` : Standard Error
* `TSTAT.[PHENO]` : T-statistic
* `LOG10P.[PHENO]` : Log10 p-value

## Full Usage 

The full usage of `qpgentools rect-assoc` can be viewed with the `--help` option:

```
$ ./qpgentools rect-assoc --help    
[bin/qpgentools rect-assoc] -- Perform rectangular association analysis for multiple phenotypes and variants

 Copyright (c) 2009-2024 by Hyun Min Kang and Adrian Tan
 Licensed under the Apache License v2.0 http://www.apache.org/licenses/

Detailed instructions of parameters are available. Ones with "[]" are in effect:

Available Options:

== Input Options ==
   --pgen            [STR: ]             : Input PLINK2 genotype file
   --psam            [STR: ]             : Input PLINK2 sample file
   --pivar           [STR: ]             : Input PLINK2 index pvar file (bgzipped and tabix)
   --pgen-list       [STR: ]             : Input file containing the list of region-specific BED files
   --pheno           [STR: ]             : Input phenotype matrix
   --sample          [STR: ]             : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
   --cov             [STR: ]             : Input covariate matrix (optional)
   --var-list        [STR: ]             : Input file containing the list of variants to be tested
   --pheno-list      [STR: ]             : Input file containing the list of phenotypes to be tested. If not provided, all phenotypes in the phenotype matrix will be tested.
   --pheno-format    [STR: regenie]      : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
   --cov-format      [STR: regenie]      : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'

== Output options ==
   --out             [STR: ]             : Output prefix

== Analysis options ==
   --rint-before-adj [FLG: OFF]          : Perform rank-based inverse normal transformation before covariate adjustment (default: false)
   --rint-after-adj  [FLG: OFF]          : Perform rank-based inverse normal transformation after covariate adjustment (default: false)
   --max-chunk-vars  [INT: 100]          : Maximum number of variants to store at once in memory (default: 100)
   --min-maf         [FLT: 1.0e-10]      : Minimum minor allele frequency to include a variant in the analysis (default: 0.0)
   --min-mac         [FLT: 1.00]         : Minimum minor allele count to include a variant in the analysis (default: 0.0)

== Auxiliary options ==
   --jump-thres-bp   [INT: 1000000]      : Jump threshold in base pairs for the variant index (default: 1000000)
   --icol-pivar-idx  [INT: 9]            : 1-based column index for the variant ID in the pvar file (default: 9)
   --icol-pheno-id   [INT: 1]            : 1-based column index for the phenotype ID in the phenotype file (default: 1)
   --icol-cov-id     [INT: 1]            : 1-based column index for the covariate ID in the phenotype file (default: 1)


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
