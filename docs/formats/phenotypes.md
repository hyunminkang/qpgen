# Phenotype File Format for `qpgentools`

## Overview

`qpgentools` accepts phenotype or covariatefiles in the following formats:

1. Regenie format
2. TensorQTL format
3. Sample-column TSV format
4. Sample-row TSV format

## Specifying phenotype or covariate files

In `qpgentools`, many commands that accepted this modified PLINK2 format, including 
`qpgentools pair-assoc`, `qpgentools pair-prs`, and `qpgentools rect-assoc`, can accept the phenotype or covariate files in the above formats. The formats can be specified using the `--pheno-format` or `--cov-format` options.

## Regenie format

The Regenie phenotype format (with `--pheno-format regenie`) is a TSV file following convention:

* The first column is `FID`, representing family ID of the sample (PLINK2 convention)
* The second column is `IID`, representing individual ID of the sample (PLINK2 convention)
* From the third column, the remaining columns are the phenotype names. 

Each row represents a sample, and each column represents a phenotype.

Here is an example:

```
FID IID pheno1 pheno2 pheno3
1 1 1.2 3.4 5.6
2 2 2.3 4.5 6.7
...
```

Please refer to the [Regenie documentation](https://rgcgithub.github.io/regenie/options/#phenotype-file-format) for more details.

## TensorQTL format

The TensorQTL phenotype format is (with `--pheno-format tensorqtl`) a TSV file that represents the phenotype in a sample-column format.

* The first column is `#chr`, representing the chromosome of the molecular trait (e.g. gene)
* The second column is `start`, representing the base position (1-based) of the molecular trait (e.g. gene)
* The third column is `end`, representing the base position (1-based) of the molecular trait (e.g. gene)
* The fourth column is `phenotype_id`, representing the ID of the molecular trait (e.g. gene)
* From the fifth column, the remaining columns are the sample IDs. 

Each row represents a molecular trait, and each column represents a sample.

Here is an example:

```
#chr    start   end     phenotype_id    sample1 sample2 sample3
1   100000  100001  gene1   1.2   3.4   5.6
2   200000  200001  gene2   2.3   4.5   6.7
...
```

## Sample-row TSV format

The sample-row TSV format (with `--pheno-format tsv-sample-row`) is a TSV file that represents the phenotype in a sample-row format. This is very similar to the Regenie format, but instead of the `FID` and `IID` columns are expected for family and individual ID, only one sample ID column is expected. The rest of the columns are the phenotype names.

Here is an example:

```
SampleID pheno1 pheno2 pheno3
sample1 1.2 3.4 5.6
sample2 2.3 4.5 6.7
...
```

## Sample-column TSV format

The sample-column TSV format (with `--pheno-format tsv-sample-column`) is a TSV file that represents the phenotype in a sample-column format. This is very similar to the TensorQTL format, but instead of the `#chr`, `start`, `end`, and `phenotype_id` columns are expected for chromosome, base position, base position, and molecular trait ID, only one phenotype ID column is expected. The rest of the columns are the sample IDs.

Here is an example:

```
TraitID sample1 sample2 sample3
trait1  1.2   3.4   5.6
trait2  2.3   4.5   6.7
...
```

 


