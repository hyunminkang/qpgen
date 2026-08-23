#ifndef __SUSIE_UTILS_H
#define __SUSIE_UTILS_H

// SuSiE ("Sum of Single Effects") fine-mapping via the IBSS algorithm, in Eigen.
//
// Core algorithm reimplemented from:
//   G. Wang, A. Sarkar, P. Carbonetto, M. Stephens (2020),
//   "A simple new approach to variable selection in regression, with application
//    to genetic fine mapping", JRSS-B 82(5):1273-1300.
// Output field names mirror susieR's susie() fit object (alpha, mu, mu2,
// lbf_variable, lbf, pip, V, sigma2, elbo, ...) plus credible sets.

#include "Eigen/Dense"
#include <vector>
#include <cstdint>

namespace susie {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

// ---- Fit object, mirroring susieR's susie() ------------------------------
struct SusieFit {
    MatrixXd alpha;        // L x p   posterior inclusion prob within each single effect
    MatrixXd mu;           // L x p   posterior mean | inclusion
    MatrixXd mu2;          // L x p   posterior 2nd moment | inclusion
    MatrixXd lbf_variable; // L x p   log Bayes factor per variable & effect
    VectorXd lbf;          // L       log Bayes factor per single effect
    VectorXd pip;          // p       marginal PIP = 1 - prod_l(1 - alpha_lj)
    VectorXd Xr;           // n       X %*% colSums(alpha*mu)
    VectorXd fitted;       // n       fitted values (= Xr; intercept = 0)
    VectorXd V;            // L       prior variance of each single effect
    double   sigma2;       // residual variance
    double   intercept;    // intercept (0 when data pre-centered)
    std::vector<double> elbo; // objective at each outer iteration
    int      niter;
    bool     converged;
    // SuSiE-inf / SuSiE-ash (unmappable effects) extras.
    //   tau2   : infinitesimal variance (inf); or sum(pi_k*sa2_k)*sigma2 (ash, for reporting)
    //   theta  : p-vector of posterior mean unmappable effects (standardized-X scale)
    //   ash_sa2: K component prior variances on the sa2 grid (ash only; empty otherwise)
    //   ash_pi : K mixture weights (ash only; sum to 1, first entry is null component)
    // For standard SuSiE and inf, ash_sa2/ash_pi are empty and theta may be empty.
    double   tau2 = 0.0;
    VectorXd theta;
    VectorXd ash_sa2;
    VectorXd ash_pi;
};

struct SusieOptions {
    int    L            = 10;
    int    max_iter     = 100;
    double tol          = 1e-3;
    bool   estimate_prior_variance    = true;
    bool   estimate_residual_variance = true;
    double scaled_prior_variance = 0.2;   // prior var = var(y)*this (init/fixed)
    double prior_v_min  = 1e-9;
    double prior_v_max  = 1e3;
    bool   standardize  = true;           // scale X columns to unit variance (a copy)
    enum ConvergenceMethod { ELBO, PIP } convergence_method = ELBO;
    // Unmappable-effects model.
    //   NONE : standard SuSiE.
    //   INF  : SuSiE-inf; theta_j ~ N(0, tau2) with MoM variance components.
    //   ASH  : SuSiE-ash; theta_j ~ sum_k pi_k * N(0, sa2_k * sigma2) with a fixed
    //          log-spaced grid (Mr.ASH-style). Simplified port that skips the
    //          LD-masking heuristics and slot-activity model to keep it fast.
    enum UnmappableEffects { NONE, INF, ASH } unmappable_effects = NONE;
    // ASH-only knobs (ignored otherwise).
    int    ash_K            = 10;        // number of prior-variance grid components (incl. null)
    double ash_sd_mult      = 2.0;       // successive sa2 components are sd_mult^2 x apart
    int    ash_inner_max    = 200;       // max coord-ascent sweeps per outer iter
    double ash_inner_tol    = 1e-3;      // loose inner tolerance; auto-tightened near outer conv.
    // Reduction-test hook: fix pi to a specific vector and skip EM updates.
    //   empty -> normal ash (EM-updated pi).
    //   length K vector -> pi held at these values throughout; used to prove
    //     ash reduces to standard SuSiE when pi=(1,0,...,0), and to SuSiE-inf
    //     when pi has all mass on a single non-null grid component.
    std::vector<double> ash_fix_pi;
    // Reduction-test hook: fix the sa2 grid and skip the data-driven rebuild.
    //   empty -> data-driven grid (see fit_susie_ash step 3).
    //   length matching ash_fix_pi -> sa2 grid held at these values.
    std::vector<double> ash_fix_sa2;
};

// A 95% (or requested coverage) credible set for a single effect.
struct CredibleSet {
    int32_t effect_index;            // 0-based single-effect index (l)
    std::vector<int32_t> variables;  // 0-based variant indices, ordered by descending alpha
    std::vector<double>  pips;        // marginal PIP of each member (parallel to variables)
    double coverage;                 // achieved cumulative alpha within the set
    double requested_coverage;
    double purity;                   // minimum absolute pairwise correlation among members
    double lbf;                      // log Bayes factor of this single effect
    double V;                        // estimated prior variance of this single effect
};

// Fit standard SuSiE. Preconditions: X columns centered (and, if opt.standardize
// is false, already scaled); y centered; no missing; no intercept column.
SusieFit fit_susie(const MatrixXd& X, const VectorXd& y, const SusieOptions& opt = SusieOptions());

// Fit SuSiE-inf (unmappable_effects = "inf"). Same preconditions as fit_susie:
// X columns centered+standardized, y centered. Adds an infinitesimal effect via
// the eigenspace (thin SVD / Gram) Omega = (tau2 XX' + sigma2 I)^{-1} formulation,
// estimating (sigma2, tau2) by method of moments and theta as its BLUP.
SusieFit fit_susie_inf(const MatrixXd& X, const VectorXd& y, const SusieOptions& opt = SusieOptions());

// Fit SuSiE-ash (simplified port of susieR unmappable_effects = "ash"). theta_j
// has a scale-mixture-of-normals prior sum_k pi_k * N(0, sa2_k * sigma2) fit by
// Mr.ASH coordinate ascent. The grid is fixed and log-spaced (first component
// is the null point mass); pi and sigma2 are estimated by EM. Skips susieR's
// LD-masking and slot-activity heuristics for speed.
SusieFit fit_susie_ash(const MatrixXd& X, const VectorXd& y, const SusieOptions& opt = SusieOptions());

// Extract credible sets from a fit. X_std must be the (standardized) design matrix
// actually used by fit_susie, so correlations/purity are computed consistently.
// Sets whose purity < min_abs_corr are dropped; duplicate sets are removed.
std::vector<CredibleSet> susie_get_cs(const SusieFit& fit, const MatrixXd& X_std,
                                      double coverage = 0.95, double min_abs_corr = 0.5);

// Bundle of everything a caller typically wants.
struct SusieResult {
    SusieFit fit;                     // includes pip, lbf, V, sigma2, elbo, niter, converged
    std::vector<CredibleSet> cs;
};

// High-level entry point requested by callers: run SuSiE fine-mapping of a single
// (covariate-adjusted, non-missing) phenotype against a genotype matrix.
//   pheno_vec : length-n phenotype (residuals after covariate adjustment); centered internally
//   geno_mat  : n x p genotype matrix (e.g. mean-centered dosages); a standardized copy is used
// Additional behavior is controlled through SusieOptions and the CS parameters.
SusieResult simple_susie_without_missing(const VectorXd& pheno_vec,
                                         const MatrixXd& geno_mat,
                                         const SusieOptions& opt = SusieOptions(),
                                         double coverage = 0.95,
                                         double min_abs_corr = 0.5);

} // namespace susie

#endif // __SUSIE_UTILS_H
