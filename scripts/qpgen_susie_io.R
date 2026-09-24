## qpgen_susie_io.R -- read the binary summary-data files written by
##   qpgentools region-assoc --out-suff / --out-rss
## and hand them to susieR (https://github.com/stephenslab/susieR).
##
## Usage:
##   source("scripts/qpgen_susie_io.R")
##   ss  <- read_qpgen_suff("out.suff.bin.gz")   # X'X, X'y, y'y, n
##   rss <- read_qpgen_rss("out.rss.bin.gz")     # z, R, n, bhat, shat, var_y
##   fit <- susie_from_qpgen(ss,  trait = "ENSG00000164308.17", L = 10)  # susie_ss(XtX, Xty, yty, n)
##   fit <- susie_from_qpgen(rss, trait = "ENSG00000164308.17", L = 10)  # susie_rss(bhat, shat, R, n, var_y)
##   fit$sets$cs; susieR::susie_get_pip(fit)
##
## Both files are gzip (BGZF) streams readable with gzfile(); no extra
## packages are needed to read them. Matrices are stored column-major in
## little-endian float64, so readBin() + matrix() reconstructs them directly.

.qpgen_read_header <- function(con, expected_magic) {
  magic <- readChar(con, 8, useBytes = TRUE)
  if (!identical(magic, expected_magic))
    stop(sprintf("unexpected magic '%s' (expected '%s')", magic, expected_magic))
  version <- readBin(con, "integer", 1, size = 4, endian = "little")
  if (version != 1L) stop("unsupported format version ", version)
  dims <- readBin(con, "integer", 4, size = 4, endian = "little")
  n <- dims[1]; p <- dims[2]; K <- dims[3]; n_cov <- dims[4]
  region   <- readBin(con, "character", 1)
  traits   <- readBin(con, "character", K)
  variants <- readBin(con, "character", p)
  pos      <- readBin(con, "integer", p, size = 4, endian = "little")
  af       <- readBin(con, "double", p, size = 8, endian = "little")
  list(version = version, n = n, p = p, K = K, n_cov = n_cov, region = region,
       traits = traits, variants = variants, pos = pos, af = af)
}

.qpgen_read_matrix <- function(con, nrow, ncol, rownames = NULL, colnames = NULL) {
  m <- matrix(readBin(con, "double", as.numeric(nrow) * as.numeric(ncol),
                      size = 8, endian = "little"), nrow = nrow, ncol = ncol)
  dimnames(m) <- list(rownames, colnames)
  m
}

## Sufficient statistics: list(XtX [p x p], Xty [p x K], yty [K], n, ...)
read_qpgen_suff <- function(path) {
  con <- gzfile(path, "rb"); on.exit(close(con))
  h <- .qpgen_read_header(con, "QPGNSUFF")
  h$XtX <- .qpgen_read_matrix(con, h$p, h$p, h$variants, h$variants)
  h$Xty <- .qpgen_read_matrix(con, h$p, h$K, h$variants, h$traits)
  h$yty <- setNames(readBin(con, "double", h$K, size = 8, endian = "little"), h$traits)
  h
}

## RSS summary statistics: list(R [p x p], z/bhat/shat [p x K], var_y [K], n, ...)
read_qpgen_rss <- function(path) {
  con <- gzfile(path, "rb"); on.exit(close(con))
  h <- .qpgen_read_header(con, "QPGN_RSS")
  h$R     <- .qpgen_read_matrix(con, h$p, h$p, h$variants, h$variants)
  h$z     <- .qpgen_read_matrix(con, h$p, h$K, h$variants, h$traits)
  h$bhat  <- .qpgen_read_matrix(con, h$p, h$K, h$variants, h$traits)
  h$shat  <- .qpgen_read_matrix(con, h$p, h$K, h$variants, h$traits)
  h$var_y <- setNames(readBin(con, "double", h$K, size = 8, endian = "little"), h$traits)
  h
}

## Fit susieR directly on the sufficient statistics of one trait:
##   susie_ss(XtX, Xty, yty, n)        in susieR >= 0.14 (exported as susie_ss)
##   susie_suff_stat(XtX, Xty, yty, n) in susieR <= 0.12 (CRAN)
## Extra arguments (L, coverage, min_abs_corr, estimate_residual_variance, ...)
## pass through. Coefficients are on the per-ALT-allele scale of the input.
susie_from_qpgen_suff <- function(ss, trait = 1, ...) {
  if (is.character(trait)) trait <- match(trait, ss$traits)
  if (is.na(trait)) stop("trait not found in file")
  Xty <- ss$Xty[, trait]; yty <- unname(ss$yty[trait])
  if (exists("susie_ss", envir = asNamespace("susieR")))
    susieR::susie_ss(XtX = ss$XtX, Xty = Xty, yty = yty, n = ss$n, ...)
  else
    susieR::susie_suff_stat(XtX = ss$XtX, Xty = Xty, yty = yty, n = ss$n, ...)
}

## Fit susieR on the RSS statistics of one trait via susie_rss(). Prefer the
## sufficient-statistics file and susie_from_qpgen_suff() when available.
## Supplying bhat/shat/var_y (rather than z alone) lets susie_rss rebuild X'X,
## X'y and y'y from R, n and var_y, but susie_rss() fixes the residual
## variance by default; pass estimate_residual_variance = TRUE to reproduce
## the sufficient-statistics fit (matches susie_ss() to ~1e-9 on in-sample R).
susie_from_qpgen_rss <- function(rss, trait = 1, ...) {
  if (is.character(trait)) trait <- match(trait, rss$traits)
  if (is.na(trait)) stop("trait not found in file")
  susieR::susie_rss(bhat = rss$bhat[, trait], shat = rss$shat[, trait], R = rss$R,
                    n = rss$n, var_y = unname(rss$var_y[trait]), ...)
}

## Dispatch on the file type: sufficient statistics -> susie_ss(); RSS -> susie_rss().
susie_from_qpgen <- function(obj, trait = 1, ...) {
  if (!is.null(obj$XtX)) susie_from_qpgen_suff(obj, trait, ...)
  else susie_from_qpgen_rss(obj, trait, ...)
}

## Derive the RSS quantities (bhat, shat, z, R, var_y) from sufficient
## statistics, using the same formulas as region-assoc's marginal test
## (centered X and y, no intercept, df = n - 2). Provided for convenience
## (e.g. to feed other RSS-based tools); not needed to run susieR.
suff_to_rss <- function(ss, trait = 1) {
  if (is.character(trait)) trait <- match(trait, ss$traits)
  d    <- diag(ss$XtX)
  Xty  <- ss$Xty[, trait]
  bhat <- Xty / d
  sse  <- pmax(ss$yty[trait] - Xty^2 / d, 0)
  shat <- sqrt(sse / (ss$n - 2) / d)
  R    <- cov2cor(ss$XtX); R <- (R + t(R)) / 2   # bit-exact symmetry
  list(bhat = bhat, shat = shat, z = bhat / shat, R = R, n = ss$n,
       var_y = unname(ss$yty[trait]) / (ss$n - 1))
}
