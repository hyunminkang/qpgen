# Genotype File Format for `qpgentools`

## Overview

`qpgentools` uses PLINK2 format genotype files as the primary input format, but it requires a modified version of `pvar` file, which is called as "indexed pvar" format, to allow fast partial access.

## Indexed PVAR format

The indexed PVAR format is a VCF-like format to store variant information, except that it has additional column to indicate the variant index as integer. The expected header is as follows (tab-delimited):

```
CHROM   POS     ID      REF     ALT     QUAL    FILTER  INFO    INDEX
```

The first 8 columns are the standard VCF columns, and the last column is the variant index as 1-based integer. The variant index is used to quickly access the variant in the genotype file.

These files are bgzipped and tabix-indexed, so that it can be accessed with standard tabix commands or libraries. In `qpgentools`, it file is called as `pivar` file, and the bgzipped files are typically named as `[prefix].var.idx.gz`, and the tabix index file is named as `[prefix].var.idx.gz.tbi`.

## Specifying genotype files

In `qpgentools`, many commands that accepted this modified PLINK2 format, including 
`qpgentools pair-assoc`, `qpgentools pair-prs`, and `qpgentools rect-assoc`, can accept the genotype files in two ways:

1. Using `--pgen`, `--psam`, and `--pivar` options, which takes the genotype, sample, and pvar files, respectively. 
    * `--pgen` : PLINK2 binary genotype file, typically named as `[prefix].pgen`
    * `--psam` : PLINK2 binary sample file, typically named as `[prefix].psam`
    * `--pivar` : Indexed pvar file, typically named as `[prefix].var.idx.gz`

2. If the genotype files are separated by chromosomes or chinks, using `--pgen-list` option, which takes a TSV file containing the list of genotype files separated by chunk. Each line is expected to have the following format

```
[CHROM] [BEG] [END] [PLINK2_PREFIX]
```
or 

```
[CHROM] [BEG] [END] [PGEN_PATH] [PSAM_PATH] [PIVAR_PATH]
```

If `[PLINK2_PREFIX]` is provided, `[PLINK2_PREFIX].pgen`, `[PLINK2_PREFIX].psam`, and `[PLINK2_PREFIX].var.idx.gz` are used as the genotype, sample, and pvar files, respectively.

For example, if the genotype files are separated by chromosomes, the `--pgen-list` option can be written as follows:

```
1   1   1000000000  /path/to/plink/prefix_chr1
2   1   1000000000  /path/to/plink/prefix_chr2
...
22  1   1000000000  /path/to/plink/prefix_chr22
23  1   1000000000  /path/to/plink/prefix_chrX
```


