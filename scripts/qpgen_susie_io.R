## qpgen_susie_io.R -- read the binary summary-data files written by
##   qpgentools region-assoc --out-suff / --out-rss
## and hand them to susieR (https://github.com/stephenslab/susieR).
##
## Usage:
##   source("scripts/qpgen_susie_io.R")
##   ss  <- read_qpgen_suff("out.suff.bin.gz")   # X'X, X'y, y'y, n
##   rss <- read_qpgen_rss("out.rss.bin.gz")     # z, R, n, bhat, shat, var_y
##   fit <- susie_from_qpgen(rss, trait = "ENSG00000164308.17", L = 10)
##   susieR::susie_get_cs(fit); susieR::susie_get_pip(fit)
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

## Convert sufficient statistics into the RSS quantities susie_rss() takes
## (bhat, shat, R, var_y). Same formulas as region-assoc's marginal test
## (centered X and y, no intercept term, df = n - 2).
suff_to_rss <- function(ss, trait = 1) {
  d    <- diag(ss$XtX)
  Xty  <- ss$Xty[, trait]
  bhat <- Xty / d
  sse  <- pmax(ss$yty[trait] - Xty^2 / d, 0)
  shat <- sqrt(sse / (ss$n - 2) / d)
  R    <- cov2cor(ss$XtX); R <- (R + t(R)) / 2   # bit-exact symmetry
  list(bhat = bhat, shat = shat, z = bhat / shat, R = R, n = ss$n,
       var_y = unname(ss$yty[trait]) / (ss$n - 1))
}

## Run susieR on one trait from either file type. Extra arguments (L,
## coverage, min_abs_corr, estimate_residual_variance, ...) pass through
## to susieR::susie_rss(). Supplying bhat/shat/var_y (rather than z alone)
## makes susie_rss equivalent to fitting on the sufficient statistics, so
## coefficients come back on the per-ALT-allele scale of the input.
susie_from_qpgen <- function(obj, trait = 1, ...) {
  if (is.character(trait)) trait <- match(trait, obj$traits)
  if (is.na(trait)) stop("trait not found in file")
  if (!is.null(obj$XtX)) {                 # sufficient statistics file
    r <- suff_to_rss(obj, trait)
    bhat <- r$bhat; shat <- r$shat; R <- r$R; var_y <- r$var_y
  } else {                                 # RSS file
    bhat <- obj$bhat[, trait]; shat <- obj$shat[, trait]; R <- obj$R
    var_y <- obj$var_y[trait]
  }
  susieR::susie_rss(bhat = bhat, shat = shat, R = R, n = obj$n,
                    var_y = unname(var_y), ...)
}
