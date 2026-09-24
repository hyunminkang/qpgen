#ifndef __SUSIE_EXPORT_H
#define __SUSIE_EXPORT_H

// Export of region-level summary data for external fine-mapping with susieR
// (https://github.com/stephenslab/susieR).
//
// Two binary, BGZF/gzip-compressed formats are written, both readable from R
// with gzfile() + readBin() (see scripts/qpgen_susie_io.R):
//
//   * sufficient statistics  ("QPGNSUFF"): X'X, X'y, y'y, n   [primary, exact]
//       -> susie_ss(XtX, Xty, yty, n) (susieR >= 0.14; susie_suff_stat() in <= 0.12)
//   * RSS summary statistics ("QPGN_RSS"): z, R (LD correlation), n, bhat, shat, var_y
//       -> susie_rss(z = z, R = R, n = n) or susie_rss(bhat, shat, R, n, var_y);
//          the latter matches the sufficient-statistics fit only with
//          estimate_residual_variance = TRUE
//
// X is the covariate-residualized, mean-centered genotype (dosage) matrix and y
// the covariate-residualized, mean-centered phenotype, i.e. exactly the data on
// which region-assoc computes the marginal association and its built-in SuSiE.
//
// File layout (all little-endian; "str" = NUL-terminated byte string;
// matrices are column-major, as R's matrix() expects):
//
//   char[8]  magic            "QPGNSUFF" or "QPGN_RSS"
//   int32    version          1
//   int32    n                number of samples
//   int32    p                number of variants
//   int32    K                number of traits
//   int32    n_cov            number of covariates residualized out (0 if none)
//   str      region           analyzed region, CHROM:BEG-END
//   str[K]   trait_ids
//   str[p]   variant_ids      CHROM:POS:REF:ALT (ALT is the effect/dosage allele)
//   int32[p] pos              1-based base position
//   double[p] af              ALT allele frequency
//   --- QPGNSUFF body ---
//   double[p*p] XtX
//   double[p*K] Xty           column k = X'y for trait k
//   double[K]   yty
//   --- QPGN_RSS body ---
//   double[p*p] R             LD correlation matrix, cov2cor(XtX)
//   double[p*K] z             marginal z-scores (= bhat/shat)
//   double[p*K] bhat          marginal effect sizes (per ALT allele)
//   double[p*K] shat          marginal standard errors
//   double[K]   var_y         sample variance of y, y'y/(n-1)

#include "Eigen/Dense"
#include <string>
#include <vector>
#include <cstdint>

namespace susie_export {

// Metadata shared by both formats.
struct RegionMeta {
    std::string region;                    // CHROM:BEG-END
    std::vector<std::string> trait_ids;    // K
    std::vector<std::string> variant_ids;  // p, CHROM:POS:REF:ALT
    std::vector<int32_t> positions;        // p
    std::vector<double> af;                // p, ALT allele frequency
    int32_t n_samples = 0;
    int32_t n_cov = 0;
};

// Write sufficient statistics. XtX is p x p, Xty is p x K, yty has length K.
// Aborts via error() on I/O failure.
void write_suff_stats(const char* path, const RegionMeta& meta,
                      const Eigen::MatrixXd& XtX, const Eigen::MatrixXd& Xty,
                      const Eigen::VectorXd& yty);

// Write RSS summary statistics. R is p x p; z, bhat, shat are p x K; var_y has length K.
void write_rss_stats(const char* path, const RegionMeta& meta,
                     const Eigen::MatrixXd& R, const Eigen::MatrixXd& z,
                     const Eigen::MatrixXd& bhat, const Eigen::MatrixXd& shat,
                     const Eigen::VectorXd& var_y);

// Convert X'X into a correlation matrix. Columns with zero variance get R_jj = 1
// and R_jk = 0 so the result stays a valid input for susie_rss.
Eigen::MatrixXd xtx_to_corr(const Eigen::MatrixXd& XtX);

} // namespace susie_export

#endif // __SUSIE_EXPORT_H
