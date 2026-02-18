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

* `--cov` : Input covariate matrix in Regenie, TensorQTL, or TSV format. See [Phenotype file formats](../formats/phenotypes.md) for more details.
* `--sample` : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
* `--pheno-format` : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
* `--cov-format` : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
* `--skip-rint` : Skip tests based on rank-based inverse normal transformation (default: false)
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

## Expected Output

The output TSV file contains the following columns, largely resembling the Regenie output format:

* `TRAIT` : Trait name
* `CHROM` : Chromosome of the variant tested
* `GENPOS` : Base position (1-based) of the variant tested
* `ID` : Variant ID in `[CHROM]:[GENPOS]:[REF]:[ALT]` format
* `ALLELE0` : Reference allele
* `ALLELE1` : Alternative allele
* `A1FREQ`: Allele frequency of alternate allele
* `N`: Sample size
* `N_RR`: Number of RR genotypes
* `N_RA`: Number of RA genotypes
* `N_AA`: Number of AA genotypes
* `INFO`: Variant quality score (typically imputation quality)
* `BETA`: Effect size for raw phenotype association
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