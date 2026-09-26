# qpgentools region-assoc

## Summary 

`qpgentools region-assoc` loads **every variant in one genomic region** at once and tests it
against one or more traits. Unlike [`rect-assoc`](rect_assoc.md), which takes an explicit variant
list and streams it in chunks, `region-assoc` keeps the whole region in memory, which makes it
the entry point for region-based methods that need the full local LD structure --- in particular
**SuSiE fine-mapping** (`--susie`).

A typical running example command is given below:

```bash
qpgentools region-assoc --pgen-list [list] --pheno [pheno] --cov [cov] \
    --region chr1:1000000-2000000 --traits ENSG00000187634 --out [out_prefix]
```

To export the region's summary data for fine-mapping in R with
[susieR](https://github.com/stephenslab/susieR) instead of (or in addition to) the built-in
fine-mapper:

```bash
qpgentools region-assoc --pgen-list [list] --pheno [pheno] --cov [cov] \
    --region chr1:1000000-2000000 --traits ENSG00000187634 --out [out_prefix] \
    --out-suff --out-rss
```

To additionally fine-map the region with SuSiE:

```bash
qpgentools region-assoc --pgen-list [list] --pheno [pheno] --cov [cov] \
    --region chr1:1000000-2000000 --traits ENSG00000187634 --out [out_prefix] \
    --susie --susie-L 10 --susie-coverage 0.95 --output-lbf
```

## Required options

* Input genotype files : Either `--pgen-list` or `--pgen`, `--psam`, `--pivar` options are required. See [Genotype file formats](../formats/genotypes.md) for more details.
* `--pheno` : Input phenotype matrix in Regenie, TensorQTL, or TSV format. See [Phenotype file formats](../formats/phenotypes.md) for more details.
* `--region` : Region to be analyzed, in `CHROM:BEG-END` format (e.g. `chr12:6534517-6538371`). The chromosome name must match the names used in the indexed pvar file.
* `--traits` or `--traitf` : Traits to be tested. `--traits` takes a comma-separated list of trait IDs; `--traitf` takes a file with one trait ID per line. Exactly one of the two must be provided.
* `--out` : Output prefix. Output file names are formed by appending the suffixes described in [Expected Output](#expected-output).

## Additional Options

### Input and sample handling

* `--cov` : Input covariate matrix in Regenie or TSV format. When provided, **both** the phenotypes and the genotypes are residualized against the covariates (see [Covariate handling](#covariate-handling)).
* `--sample` : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files.
* `--loco` : REGENIE step-1 leave-one-chromosome-out (LOCO) predictions, either a single `.loco` file (only when one trait is tested) or a REGENIE `*_pred.list` file with one `TRAIT PATH` line per trait. The prediction row for the chromosome of `--region` is subtracted from the phenotype (see [LOCO adjustment](#loco-adjustment)).
* `--pheno-format` : Format of the phenotype file (default: `regenie`). Options: `regenie`, `tensorqtl`, `tsv-sample-col`, `tsv-sample-row`.
* `--cov-format` : Format of the covariate file (default: `regenie`). Options: `regenie`, `tsv-sample-col`, `tsv-sample-row`.
* `--colname-pheno-sample` / `--colname-geno-sample` : When `--sample` is provided, the column names holding the phenotype-side and genotype-side sample IDs (defaults: `pheno` and `geno`).
* `--icol-pivar-idx` : 1-based column index for the variant index in the indexed pvar file (default: 9).
* `--icol-pheno-id` / `--icol-cov-id` : 1-based column index for the phenotype/covariate ID (default: 1).
* `--offset-pheno` / `--offset-cov` : The number of leading index columns (i.e. not containing values) in the phenotype/covariate file (default: 1).

### Variant filters

* `--min-af` / `--max-af` : Minimum / maximum ALT allele frequency for a variant to be included (defaults: 0.0 and 1.0).
* `--min-ac` / `--max-ac` : Minimum / maximum ALT allele count for a variant to be included (defaults: 0.0 and 1e9).
* `--max-allowed-vars` : Maximum number of variants that may be held in memory for the region (default: 1000000). The run aborts if the region contains at least this many variants, since the whole region must fit in memory.
* `--jump-thres-bp` : Jump threshold in base pairs for the variant index (default: 1000000).

Monomorphic variants (AC = 0 or AC = AN) in the analyzed samples are always dropped, regardless of the filters above.

### Transformations

* `--rint-before-adj` : Perform rank-based inverse normal transformation before covariate adjustment (default: false).
* `--rint-after-adj` : Perform rank-based inverse normal transformation after covariate adjustment (default: false).

RINT uses Blom offsets, \(\Phi^{-1}\big((r - 3/8)/(n + 1/4)\big)\) with average ranks \(r\) for ties, as REGENIE's `--apply-rint` and `--apply-rerint` do.

### Test statistic

* `--score` : Use the REGENIE-style score test instead of the default Wald t-test (default: false). The association output then has `ZSTAT` columns instead of `TSTAT`.

The two tests give the same `BETA = x'y / x'x` but differ in how the residual variance is
estimated:

* **Wald t-test (default).** The residual variance is estimated after fitting the variant,
  `SE = sqrt((y'y - (x'y)^2 / x'x) / (N - 2) / x'x)`, and `LOG10P` comes from a t-distribution
  with `N - 2` degrees of freedom.
* **Score test (`--score`).** The residual variance is estimated under the null model,
  `sigma0^2 = y'y / (N - n_cov)`, where `n_cov` counts the intercept and the covariates. Then
  `SE = sqrt(sigma0^2 / x'x)`, `ZSTAT = BETA / SE`, and `LOG10P` comes from the standard normal.
  This matches REGENIE step 2 for quantitative traits.

For a variant that explains a fraction `R^2` of the phenotype variance, the Wald `SE` is smaller
than the score `SE` by about `sqrt(1 - R^2)`. The difference is negligible for most variants but
visible at strong hits.

### Output file names

* `--assoc-suffix` : Suffix for the marginal association output (default: `.assoc.tsv.gz`).
* `--susie-cs-suffix` : Suffix for the SuSiE credible-set output (default: `.susie.cs.tsv.gz`).
* `--susie-lbf-suffix` : Suffix for the SuSiE per-variant log Bayes factor output (default: `.susie.lbf.tsv.gz`).
* `--suff-suffix` : Suffix for the sufficient-statistics export written by `--out-suff` (default: `.suff.bin.gz`).
* `--rss-suffix` : Suffix for the RSS summary-statistics export written by `--out-rss` (default: `.rss.bin.gz`).

Suffixes ending in `.gz` produce bgzipped output; any other suffix produces plain text. The two
binary exports are always bgzipped regardless of suffix.

### Summary-data export for susieR

* `--out-suff` : Write the **sufficient statistics** `X'X`, `X'y`, `y'y` and `n` to `[out_prefix].suff.bin.gz` (default: off).
* `--out-rss` : Write the **RSS summary statistics** `z`, LD correlation matrix `R`, `n`, `bhat`, `shat` and `var_y` to `[out_prefix].rss.bin.gz` (default: off).

Both are computed from the same covariate-residualized, mean-centered genotype matrix `X` and
phenotype matrix `Y` used for the marginal association and for `--susie`, so fitting them with
`susieR` reproduces the built-in fine-mapping. See
the corresponding section under [Expected Output](#expected-output) for the file layout and how
to load them in R. Both files contain a `p x p` matrix of doubles, so their size grows as
`8 * p^2` bytes before compression: about 80 MB for 3,000 variants and 8 GB for 30,000. Restrict
the region (or the allele-frequency filters) accordingly.

## SuSiE fine-mapping

`region-assoc` implements SuSiE ("Sum of Single Effects") fine-mapping via the IBSS algorithm,
reimplemented in Eigen from Wang, Sarkar, Carbonetto & Stephens (2020), *JRSS-B* 82(5):1273-1300.
Field names in the output mirror `susieR`'s `susie()` fit object (`alpha`, `mu`, `mu2`,
`lbf_variable`, `lbf`, `pip`, `V`, `sigma2`), so results are directly comparable with `susieR`.

Fine-mapping runs on exactly the same data as the marginal association reported in the
`.assoc.tsv.gz` file: the same genotype matrix, the same covariate-adjusted phenotypes, and the
same set of samples. This is why the credible-set output can carry the marginal `BETA`, `SE`, and
`LOG10P` alongside the posterior quantities.

### SuSiE options

* `--susie` : Run SuSiE fine-mapping for each tested trait in the region (default: off). Without this flag only the marginal association file is written.
* `--susie-L` : Maximum number of causal single effects (default: 10). Effectively capped at the number of variants in the region.
* `--susie-max-iter` : Maximum number of IBSS iterations (default: 100).
* `--susie-tol` : Convergence tolerance for the SuSiE objective (default: 1e-3). For `--unmappable-effects inf` and `ash`, which converge on PIPs rather than the ELBO, a default of 1e-4 is used instead, matching `susieR`.
* `--susie-coverage` : Target coverage of the credible sets (default: 0.95).
* `--susie-min-abs-corr` : Minimum purity, i.e. the minimum absolute pairwise correlation among the members of a credible set, required to report it (default: 0.5). Sets below this threshold are dropped, as are duplicate sets.
* `--susie-no-standardize` : Do not standardize genotype columns to unit variance before fitting. Columns are still mean-centered.
* `--output-lbf` : Also write the per-variant log Bayes factors (one column per single effect) to `[out_prefix].susie.lbf.tsv.gz`.

### Unmappable-effects models

`--unmappable-effects` selects how effects that no single variant can capture (e.g. untyped
causal variants, or polygenic background within the region) are modeled. It matches the
`--method` argument of `run_susie_v1.r`.

| Value | Model |
|---|---|
| `none` (default) | Standard SuSiE. |
| `inf` | **SuSiE-inf**: adds an infinitesimal effect, `theta_j ~ N(0, tau2)`, with `(sigma2, tau2)` estimated by method of moments and `theta` obtained as its BLUP. Convergence is assessed on PIPs. |
| `ash` | **SuSiE-ash**: `theta_j` gets a scale-mixture-of-normals prior, `sum_k pi_k * N(0, sa2_k * sigma2)`, over a fixed log-spaced grid whose first component is the null point mass, fit by Mr.ASH-style coordinate ascent with `pi` and `sigma2` estimated by EM. This is a simplified port that skips `susieR`'s LD-masking and slot-activity heuristics for speed. Convergence is assessed on PIPs. |

For `inf` and `ash`, `tau2` is reported in the log and a per-variant `theta` column is added to
the LBF output.

`inf` and `ash` start from a thin eigendecomposition of the standardized genotype matrix,
taken on the smaller of `X'X` (p x p) or `XX'` (n x n), so the cost is
`O(n p min(n,p) + min(n,p)^3)` like `susieR`'s `svd(X)`. When CMake finds a system LAPACK/BLAS
(Accelerate on macOS; OpenBLAS/MKL on Linux, selectable with `-DBLA_VENDOR=...`), this step uses
the multithreaded library. Otherwise it falls back to Eigen's single-threaded routines, which
are several times slower for regions with thousands of variants.

* `--threads N` : Number of BLAS/LAPACK threads (default `0` = the library's own default, usually
  all cores). Supported for Accelerate, OpenBLAS, MKL, BLIS and FlexiBLAS. On a cluster, match it
  to the cores you requested, e.g. `--threads $SLURM_CPUS_PER_TASK`.

The backend actually loaded is detected at run time. It is shown in `region-assoc --help` (under
`--threads`) and logged as `SuSiE linear algebra: <backend> (threads: N)` when `inf`/`ash` runs.

As in `susieR`, single effects whose estimated prior variance is ~0 (`V <= 1e-9`) are treated
as inactive: they contribute neither to the reported PIPs nor to credible sets.

Two options exist purely for **reduction tests**, i.e. verifying that `ash` collapses onto the
simpler models, and are ignored unless `--unmappable-effects ash` is given:

* `--ash-fix-pi` : Comma-separated mixture weights held fixed (EM skipped). Setting `pi = 1,0,...,0` should reproduce `--unmappable-effects none`; putting all mass on a single non-null component should reproduce `inf`.
* `--ash-fix-sa2` : Comma-separated prior-variance grid held fixed. Must have the same length as `--ash-fix-pi`.

### Variants excluded from fine-mapping

Because a single effect shares one residual and one softmax across all variants, one bad column
would contaminate the whole fit rather than just its own coefficient. `region-assoc` therefore
zeroes out (and reports in the log) any variant column that is non-finite, monomorphic in the
analyzed samples, or fully explained by the covariates after residualization. Such variants get
zero prior weight, so their `alpha` and `pip` are 0, and the remaining variants are fit normally.
They still appear in the marginal association output.

## Covariate handling

When `--cov` is supplied, the phenotype matrix is adjusted for the covariates by linear
regression, and the genotype matrix is residualized against the *same* covariates
(Frisch-Waugh-Lovell). Residualizing only the phenotype would leave genotype variance that is
collinear with the covariates (e.g. genotype PCs) in the design, inflating `x'x`, attenuating the
correlation with the phenotype, and thus shrinking `BETA` and the test statistic. Residualizing
both sides makes the reported effect the *partial* effect of the variant, and matches the usual
`susieR` workflow where both `X` and `y` are residualized before fitting.

!!! note
    Without `--score`, the reported `LOG10P` uses `df = N - 2` and does not subtract the number of covariates. With
    a large sample size relative to the number of covariates the difference is negligible, but it
    is worth keeping in mind for small `N`.

!!! warning
    Covariate adjustment and the rank-based inverse normal transformations are currently supported
    only for phenotype and covariate matrices **without missing values**; the run aborts otherwise.

## LOCO adjustment

With `--loco`, the phenotype is adjusted for the polygenic background estimated by REGENIE
step 1, following REGENIE step 2 for quantitative traits:

1. RINT the raw phenotype if `--rint-before-adj` is set.
2. Residualize the phenotype against the covariates (or just center it without `--cov`).
3. Scale the residual to unit SD, with `SD = ||r|| / sqrt(N - n_cov)`, where `n_cov` counts the
   intercept (REGENIE's `scale_Y`).
4. Subtract the LOCO prediction for the chromosome of `--region`.
5. RINT the result if `--rint-after-adj` is set. Otherwise, multiply it back by the SD from
   step 3.

Because of step 5, `BETA`, `SE`, the exported `--out-suff`/`--out-rss` statistics, and the SuSiE
effect sizes stay in units of the covariate-adjusted phenotype, as REGENIE reports them.
`TSTAT` and `LOG10P` do not depend on this rescaling. With `--rint-after-adj`, effects are on
the RINT scale, as with REGENIE's `--apply-rerint`.

A `.loco` file has a header `FID_IID <FID>_<IID> ...` and one row per chromosome (`1`-`22`,
and `23` for X). The region's chromosome is matched after dropping a `chr` prefix, with X
mapped to `23`. Sample IDs are matched to the `FID_IID` columns exactly, then as `<ID>_<ID>`,
and then as a unique `<FID>_<ID>` suffix. Samples without a non-missing LOCO prediction for
every tested trait are excluded, and a warning reports how many. In a `*_pred.list` file,
relative paths that do not exist as given are resolved against the list file's directory.

```bash
qpgentools region-assoc --pgen-list [list] --pheno [pheno] --cov [cov] \
    --traits T1,T2 --region chr2:1000000-2000000 \
    --loco [regenie_step1]_pred.list --susie --out [out_prefix]
```

!!! warning
    `--loco` currently requires a phenotype matrix without missing values.

## Expected Output

### `[out_prefix].assoc.tsv.gz` --- marginal association

One row per variant in the region that passed the filters, in a wide format where the statistics
for each trait are appended horizontally.

* `#CHROM` : Chromosome
* `GENPOS` : Base position (1-based)
* `ID` : Variant ID in `[CHROM]:[GENPOS]:[REF]:[ALT]` format
* `ALLELE0` : Reference allele
* `ALLELE1` : Alternative allele
* `A1FREQ` : Frequency of the alternative allele (AC / AN)
* `N` : Number of samples with a called genotype (`N_RR + N_RA + N_AA`)
* `N_RR` : Count of Reference/Reference genotypes
* `N_RA` : Count of Reference/Alternative genotypes
* `N_AA` : Count of Alternative/Alternative genotypes
* `TEST` : Test type (always `Linear`)

Followed, for each trait `[TRAIT]`, by:

* `BETA.[TRAIT]` : Effect size
* `SE.[TRAIT]` : Standard error
* `TSTAT.[TRAIT]` : Wald t-statistic, or `ZSTAT.[TRAIT]` : score z-statistic with `--score`
* `LOG10P.[TRAIT]` : -log10 p-value of the selected test

### `[out_prefix].susie.cs.tsv.gz` --- credible sets

Written only with `--susie`. One row per (trait, credible set, variant) triple, i.e. only variants
that belong to a reported credible set appear here.

* `#trait` : Trait ID
* `variant` : Variant ID in `[CHROM]:[POS]:[REF]:[ALT]` format
* `pip` : Marginal posterior inclusion probability, `1 - prod_l (1 - alpha_lj)`
* `cs_id` : 1-based credible set index within this trait
* `alpha` : Posterior inclusion probability of this variant *within this single effect*
* `region` : Analyzed region, formatted as `CHROM_BEG_END`
* `n_region_vars` : Number of variants analyzed in this region (identical on every row of the file)
* `cs_size` : Number of variants in this credible set
* `cs_lbf` : Log Bayes factor (natural log) of the single effect defining this credible set, i.e. the prior-weighted average BF over all variants in the region. Constant across all rows of a credible set
* `var_lbf` : Log Bayes factor (natural log) of *this variant* under that single effect. Related to `alpha` by `var_lbf - cs_lbf = log(alpha) - log(prior)`
* `mu` : Posterior mean effect, conditional on inclusion
* `mu2` : Posterior second moment of the effect, conditional on inclusion
* `af` : Alternative allele frequency
* `n` : Sample size used in the marginal test
* `beta` : Marginal effect size (matches `BETA.[TRAIT]` in the association file)
* `se` : Marginal standard error (matches `SE.[TRAIT]`, so it follows `--score`)
* `log10p` : Marginal -log10 p-value (matches `LOG10P.[TRAIT]`, so it follows `--score`)

A trait with no credible set surviving the coverage and purity thresholds contributes no rows.

### `[out_prefix].susie.lbf.tsv.gz` --- per-variant log Bayes factors

Written only with `--susie --output-lbf`. One row per (trait, variant) pair for **every** variant
in the region, whether or not it belongs to a credible set.

* `#trait` : Trait ID
* `variant` : Variant ID
* `region` : Analyzed region, formatted as `CHROM_BEG_END`
* `n_region_vars` : Number of variants analyzed in this region (identical on every row of the file)
* `af` : Alternative allele frequency
* `pip` : Marginal posterior inclusion probability
* `theta` : Posterior mean unmappable effect on the standardized-genotype scale (equivalent to `susieR`'s `fit$theta`). **Present only with `--unmappable-effects inf` or `ash`.**
* `lbf.L1` ... `lbf.L[L]` : Log Bayes factor of this variant under each single effect, where `L` is `min(--susie-L, number of variants)`

### `[out_prefix].suff.bin.gz` and `[out_prefix].rss.bin.gz` --- summary-data export for susieR

Written only with `--out-suff` / `--out-rss`. These are gzip (BGZF) compressed binary files
that need no extra R packages to read: `gzfile()` + `readBin()` is enough. A ready-made reader is
provided in [`scripts/qpgen_susie_io.R`](https://github.com/hyunminkang/qpgen/blob/main/scripts/qpgen_susie_io.R).

**Quantities.** `X` is the `n x p` matrix of ALT-allele dosages, mean-imputed for missing
genotypes, mean-centered, and residualized against the covariates when `--cov` is given. `y` is
the covariate-residualized, mean-centered phenotype (one column per trait, `K` traits).

| File | Contents | susieR entry point |
|---|---|---|
| `.suff.bin.gz` | `XtX = X'X` (`p x p`), `Xty = X'y` (`p x K`), `yty = y'y` (`K`), `n` | `susie_ss(XtX, Xty, yty, n)` in susieR >= 0.14 (`susie_suff_stat()` in CRAN susieR <= 0.12) |
| `.rss.bin.gz` | `R = cov2cor(X'X)` (`p x p`), `z = bhat/shat`, `bhat`, `shat` (`p x K`), `var_y = y'y/(n-1)` (`K`), `n` | `susie_rss(z = z, R = R, n = n)` or `susie_rss(bhat, shat, R, n, var_y)` |

The sufficient-statistics file is the primary export. `X'X`, `X'y`, `y'y` and `n` are the exact
individual-level likelihood, and `susie_ss()` on them reproduces the built-in `--susie` fit
(same `sigma2`, credible sets and PIPs). Use this file whenever you fit susieR yourself.

The RSS file exists for tools that take marginal statistics plus an LD matrix. `bhat`, `shat` and
`z` in it are the same values as `BETA`, `SE` and `TSTAT` (or `ZSTAT` with `--score`) in the association file, and `R` is the
in-sample LD of the same residualized data. Note that `susie_rss()` is **not** equivalent to the
sufficient-statistics fit under its defaults: it fixes the residual variance
(`estimate_residual_variance = FALSE`), which on the test data changed `sigma2` from 0.937 to
1.285 and moved PIPs by up to 0.05. With `estimate_residual_variance = TRUE` and `bhat`, `shat`,
`var_y` and `n` supplied, `susie_rss()` rebuilds `X'X`, `X'y` and `y'y` internally and matches
`susie_ss()` to about 1e-9. Passing `z` alone fixes `sigma2 = 1` and gives the RSS-model fit.

**Binary layout** (little-endian; strings are NUL-terminated; matrices are column-major, as R's
`matrix()` expects):

```
char[8]   magic       "QPGNSUFF" (suff) or "QPGN_RSS" (rss)
int32     version     1
int32     n           number of samples
int32     p           number of variants
int32     K           number of traits
int32     n_cov       number of covariates residualized out (0 without --cov)
str       region      CHROM:BEG-END
str[K]    trait_ids
str[p]    variant_ids CHROM:POS:REF:ALT (ALT is the dosage allele)
int32[p]  pos         1-based position
double[p] af          ALT allele frequency
--- suff body ---              --- rss body ---
double[p*p] XtX                double[p*p] R
double[p*K] Xty                double[p*K] z
double[K]   yty                double[p*K] bhat
                               double[p*K] shat
                               double[K]   var_y
```

**Loading in R.** With the helper script:

```r
library(susieR)
source("scripts/qpgen_susie_io.R")

ss  <- read_qpgen_suff("out.suff.bin.gz")   # list(XtX, Xty, yty, n, variants, traits, af, ...)
fit <- susie_from_qpgen(ss, trait = "ENSG00000187634", L = 10)   # = susie_ss(XtX, Xty, yty, n, L = 10)
fit$sets$cs                                  # credible sets (indices into ss$variants)
head(sort(susie_get_pip(fit), decreasing = TRUE))

rss  <- read_qpgen_rss("out.rss.bin.gz")    # list(R, z, bhat, shat, var_y, n, ...)
fit2 <- susie_from_qpgen(rss, trait = 1, L = 10,   # = susie_rss(bhat, shat, R, n, var_y, ...)
                         estimate_residual_variance = TRUE)   # needed to match fit
```

`susie_from_qpgen()` dispatches on the file type: `susie_ss()` for the sufficient statistics
(`susie_suff_stat()` on old susieR), `susie_rss()` for the RSS file. Any further arguments
(`L`, `coverage`, `min_abs_corr`, `estimate_residual_variance`, ...) pass through. Use
`fit$sets` for credible sets, or `susie_get_cs(fit, Xcorr = cov2cor(ss$XtX))`; calling
`susie_get_cs()` without `Xcorr` skips the purity filter.

Without the helper, the RSS file reads as:

```r
con <- gzfile("out.rss.bin.gz", "rb")
stopifnot(readChar(con, 8, useBytes = TRUE) == "QPGN_RSS")
version <- readBin(con, "integer", 1, size = 4, endian = "little")
d <- readBin(con, "integer", 4, size = 4, endian = "little"); n <- d[1]; p <- d[2]; K <- d[3]
region   <- readBin(con, "character", 1)
traits   <- readBin(con, "character", K)
variants <- readBin(con, "character", p)
pos      <- readBin(con, "integer", p, size = 4, endian = "little")
af       <- readBin(con, "double", p)
R    <- matrix(readBin(con, "double", p * p), p, p, dimnames = list(variants, variants))
z    <- matrix(readBin(con, "double", p * K), p, K, dimnames = list(variants, traits))
bhat <- matrix(readBin(con, "double", p * K), p, K)
shat <- matrix(readBin(con, "double", p * K), p, K)
var_y <- readBin(con, "double", K)
close(con)
fit <- susie_rss(bhat = bhat[, 1], shat = shat[, 1], R = R, n = n, var_y = var_y[1], L = 10)
```

The `.suff.bin.gz` file has the same header followed by `XtX` (`p*p` doubles), `Xty` (`p*K`)
and `yty` (`K`), to be passed as `susie_ss(XtX, Xty[, 1], yty[1], n, L = 10)`. Python users can read the same layout with `gzip.open()` and
`numpy.frombuffer(..., dtype="<f8").reshape((p, p), order="F")`.

!!! note
    The LD matrix `R` is the in-sample correlation of the *covariate-residualized* dosages, not
    of the raw genotypes. This is what makes `susie_rss` on these files equivalent to the
    individual-level fit; it is not a drop-in reference LD panel for other GWAS.

## Full Usage 

The full usage of `qpgentools region-assoc` can be viewed with the `--help` option:

```
$ ./qpgentools region-assoc --help
[bin/qpgentools region-assoc] -- Perform association analysis for a specific region and multiple phenotypes

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
   --sample               [STR: ]             : Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files
   --cov                  [STR: ]             : Input covariate matrix (optional)
   --loco                 [STR: ]             : REGENIE step-1 LOCO predictions: a .loco file (single trait) or a *_pred.list file (TRAIT PATH per line). The covariate-adjusted phenotype is scaled to unit SD, the prediction for the region's chromosome is subtracted, and (unless --rint-after-adj) the result is rescaled to the original SD so BETA/SE stay in phenotype units
   --pheno-format         [STR: regenie]      : Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'
   --cov-format           [STR: regenie]      : Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'
   --traits               [STR: ]             : Trait IDs (comma-separated) to be tested (required)
   --traitf               [STR: ]             : Input file containing trait IDs to be tested (one per line)
   --region               [STR: ]             : Region string in the format CHROM:BEG-END (required)
   --min-af               [FLT: 0.00]         : Minimum allele frequency for variants to be tested (default: 0.0)
   --max-af               [FLT: 1.00]         : Maximum allele frequency for variants to be tested (default: 1.0)
   --min-ac               [FLT: 0.00]         : Minimum allele count for variants to be tested (default: 0.0)
   --max-ac               [FLT: 1000000000.00] : Maximum allele count for variants to be tested (default: 1e9)
   --threads              [INT: 0]            : Number of threads for the BLAS/LAPACK linear algebra of SuSiE-inf/ash (default: 0 = library default). Linked backend: Apple Accelerate

== Output options ==
   --out                  [STR: ]             : Output prefix
   --assoc-suffix         [STR: .assoc.tsv.gz] : Suffix for the association output file (default: '.assoc.tsv.gz')
   --susie-cs-suffix      [STR: .susie.cs.tsv.gz] : Suffix for the SuSiE credible set output file (default: '.susie.cs.tsv.gz')
   --susie-lbf-suffix     [STR: .susie.lbf.tsv.gz] : Suffix for the SuSiE log Bayes factor output file (default: '.susie.lbf.tsv.gz')
   --out-suff             [FLG: OFF]          : Write sufficient statistics (X'X, X'y, y'y, n) for susieR as a gzipped binary file [out]<suff-suffix>
   --out-rss              [FLG: OFF]          : Write RSS summary statistics (z, LD matrix R, n, bhat, shat, var_y) for susieR::susie_rss as a gzipped binary file [out]<rss-suffix>
   --suff-suffix          [STR: .suff.bin.gz] : Suffix for the sufficient statistics output file (default: '.suff.bin.gz')
   --rss-suffix           [STR: .rss.bin.gz]  : Suffix for the RSS summary statistics output file (default: '.rss.bin.gz')

== Auxiliary options ==
   --jump-thres-bp        [INT: 1000000]      : Jump threshold in base pairs for the variant index (default: 1000000)
   --max-allowed-vars     [INT: 1000000]      : Maximum number of allowed variants to store at once in memory (default: 1000000)
   --offset-pheno         [INT: 1]            : The number of index columns (i.e. not containing values) in the phenotype files (default: 1)
   --offset-cov           [INT: 1]            : The number of index columns (i.e. not containing values) in the covariate columns (default: 1)
   --icol-pivar-idx       [INT: 9]            : 1-based column index for the variant ID in the pvar file (default: 9)
   --icol-pheno-id        [INT: 1]            : 1-based column index for the phenotype ID in the phenotype file (default: 1)
   --icol-cov-id          [INT: 1]            : 1-based column index for the covariate ID in the phenotype file (default: 1)
   --colname-pheno-sample [STR: pheno]        : When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')
   --colname-geno-sample  [STR: geno]         : When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'geno')
   --rint-before-adj      [FLG: OFF]          : Perform rank-based inverse normal transformation before covariate adjustment (default: false)
   --rint-after-adj       [FLG: OFF]          : Perform rank-based inverse normal transformation after covariate adjustment (default: false)
   --score                [FLG: OFF]          : Use the REGENIE-style score test (residual variance under the null with N - n_cov df, normal p-value) and report ZSTAT instead of the Wald TSTAT

== SuSiE fine-mapping options ==
   --susie                [FLG: OFF]          : Run SuSiE fine-mapping for each tested phenotype in the region
   --susie-L              [INT: 10]           : Maximum number of causal single effects (default: 10)
   --susie-max-iter       [INT: 100]          : Maximum number of IBSS iterations (default: 100)
   --susie-coverage       [FLT: 0.95]         : Target coverage of credible sets (default: 0.95)
   --susie-min-abs-corr   [FLT: 0.50]         : Minimum absolute correlation (purity) required to report a credible set (default: 0.5)
   --susie-tol            [FLT: 1.0e-03]      : Convergence tolerance for the SuSiE objective (default: 1e-3)
   --susie-no-standardize [FLG: OFF]          : Do not standardize genotype columns to unit variance before SuSiE
   --unmappable-effects   [STR: none]         : Unmappable-effects model for SuSiE: 'none' (standard), 'inf' (SuSiE-inf, adds an infinitesimal effect), or 'ash' (SuSiE-ash, scale-mixture prior). Matches run_susie_v1.r --method (default: none)
   --ash-fix-pi           [STR: ]             : Reduction test: comma-separated pi vector to hold ash mixture weights fixed (skips EM). Length K. Used with --unmappable-effects ash to prove ash reduces to none (pi=1,0,...,0) or inf (pi=0,...,0,1).
   --ash-fix-sa2          [STR: ]             : Reduction test: comma-separated sa2 grid to hold ash prior-variance grid fixed. Length must match --ash-fix-pi.
   --output-lbf           [FLG: OFF]          : Also write per-variant log Bayes factors (one column per single effect) when running SuSiE


NOTES:
When --help was included in the argument. The program prints the help message but do not actually run
```
