#include "susie_utils.h"
#include "qgenlib/qgen_error.h"
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <random>
#include <string>

namespace susie {

// ---------------------------------------------------------------------------
// Single-effect regression (SER) and prior-variance optimization
// ---------------------------------------------------------------------------

struct SER {
    VectorXd alpha;   // p
    VectorXd mu;      // p
    VectorXd mu2;     // p
    VectorXd lbf_var; // p
    double   lbf;     // scalar log Bayes factor for the single effect
    double   V;
};

static double ser_loglik(double v, const ArrayXd& betahat2,
                         const ArrayXd& s2, const ArrayXd& logpi) {
    ArrayXd lbf = 0.5 * (s2 / (s2 + v)).log()
                + 0.5 * (betahat2 / s2) * (v / (s2 + v));
    ArrayXd w = lbf + logpi;
    double m = w.maxCoeff();
    return m + std::log((w - m).exp().sum());
}

static double optimize_prior_variance(const ArrayXd& betahat2,
                                      const ArrayXd& s2,
                                      const ArrayXd& logpi,
                                      double lo, double hi) {
    // golden-section search on log(v), then compare against the null (v = 0)
    const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
    double a = std::log(lo), b = std::log(hi);
    double c = b - gr * (b - a), d = a + gr * (b - a);
    double fc = ser_loglik(std::exp(c), betahat2, s2, logpi);
    double fd = ser_loglik(std::exp(d), betahat2, s2, logpi);
    for (int i = 0; i < 100 && (b - a) > 1e-6; ++i) {
        if (fc > fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = ser_loglik(std::exp(c), betahat2, s2, logpi); }
        else         { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = ser_loglik(std::exp(d), betahat2, s2, logpi); }
    }
    double v_opt   = std::exp((a + b) / 2.0);
    double ll_opt  = ser_loglik(v_opt, betahat2, s2, logpi);
    double ll_null = ser_loglik(0.0,   betahat2, s2, logpi);
    return (ll_null > ll_opt) ? 0.0 : v_opt;
}

static SER single_effect_regression(const VectorXd& r, const MatrixXd& X,
                                     const ArrayXd& xtx, double sigma2,
                                     double V_init, const ArrayXd& logpi,
                                     bool estimate_v, double v_lo, double v_hi) {
    const int p = X.cols();
    ArrayXd xtr = (X.transpose() * r).array();
    ArrayXd betahat = xtr / xtx;
    ArrayXd s2 = sigma2 / xtx;
    ArrayXd betahat2 = betahat.square();

    double v = estimate_v ? optimize_prior_variance(betahat2, s2, logpi, v_lo, v_hi)
                          : V_init;

    SER out; out.V = v;
    if (v <= 0.0) {
        out.alpha   = VectorXd::Constant(p, 1.0 / p);
        out.mu      = VectorXd::Zero(p);
        out.mu2     = VectorXd::Zero(p);
        out.lbf_var = VectorXd::Zero(p);
        out.lbf     = 0.0;
        return out;
    }

    ArrayXd lbf = 0.5 * (s2 / (s2 + v)).log()
                + 0.5 * (betahat2 / s2) * (v / (s2 + v));
    ArrayXd w = lbf + logpi;
    double m = w.maxCoeff();
    ArrayXd a = (w - m).exp();
    a /= a.sum();

    double lbf_model = m + std::log((w - m).exp().sum()); // log average BF under prior

    ArrayXd post_var = 1.0 / (1.0 / s2 + 1.0 / v);
    ArrayXd mu  = (post_var / s2) * betahat;
    ArrayXd mu2 = mu.square() + post_var;

    out.alpha   = a.matrix();
    out.mu      = mu.matrix();
    out.mu2     = mu2.matrix();
    out.lbf_var = lbf.matrix();
    out.lbf     = lbf_model;
    return out;
}

// ---------------------------------------------------------------------------
// IBSS fit
// ---------------------------------------------------------------------------

SusieFit fit_susie(const MatrixXd& X, const VectorXd& y, const SusieOptions& opt) {
    const int n = X.rows();
    const int p = X.cols();
    const int L = std::max(1, std::min(opt.L, p));

    ArrayXd xtx = X.array().square().colwise().sum();
    // Columns with no variance (monomorphic, or fully explained by the covariates
    // after the Frisch-Waugh-Lovell residualization) carry no information. They must
    // be given zero prior weight rather than merely a guarded x'x: dividing by a
    // sentinel x'x yields betahat2/s2 = 0/0 = NaN, and because every variable shares
    // the same softmax and residual, a single NaN propagates to the whole fit
    // (all-NaN alpha/lbf, sigma2 frozen at its initial value, no credible sets).
    int n_active = 0;
    for (int j = 0; j < p; ++j) {
        if (std::isfinite(xtx(j)) && xtx(j) > 0.0) ++n_active;
    }
    if (n_active == 0) {
        error("fit_susie(): all %d genotype columns have zero variance; nothing to fine-map", p);
    }
    ArrayXd logpi(p);
    const double logpi_active = -std::log((double)n_active);
    for (int j = 0; j < p; ++j) {
        if (std::isfinite(xtx(j)) && xtx(j) > 0.0) {
            logpi(j) = logpi_active;
        } else {
            xtx(j)   = 1.0;                                        // keep the arithmetic finite
            logpi(j) = -std::numeric_limits<double>::infinity();    // prior weight 0 -> alpha 0
        }
    }
    double var_y = (y.array() - y.mean()).square().sum() / std::max(1, n);

    MatrixXd alpha = MatrixXd::Zero(L, p);
    MatrixXd mu    = MatrixXd::Zero(L, p);
    MatrixXd mu2   = MatrixXd::Zero(L, p);
    MatrixXd lbf_variable = MatrixXd::Zero(L, p);
    VectorXd lbf   = VectorXd::Zero(L);
    VectorXd V     = VectorXd::Constant(L, opt.scaled_prior_variance * var_y);
    double sigma2  = var_y > 0.0 ? var_y : 1.0;

    VectorXd b_bar = VectorXd::Zero(p);

    SusieFit fit; fit.converged = false; fit.niter = 0;
    double obj_prev = -std::numeric_limits<double>::infinity();
    MatrixXd alpha_prev = alpha;

    for (int iter = 0; iter < opt.max_iter; ++iter) {
        for (int l = 0; l < L; ++l) {
            VectorXd b_l = (alpha.row(l).array() * mu.row(l).array()).matrix();
            VectorXd r = y - X * (b_bar - b_l);   // residual excluding effect l
            SER ser = single_effect_regression(r, X, xtx, sigma2, V(l), logpi,
                                               opt.estimate_prior_variance,
                                               opt.prior_v_min, opt.prior_v_max);
            alpha.row(l)        = ser.alpha.transpose();
            mu.row(l)           = ser.mu.transpose();
            mu2.row(l)          = ser.mu2.transpose();
            lbf_variable.row(l) = ser.lbf_var.transpose();
            lbf(l)              = ser.lbf;
            V(l)                = ser.V;
            VectorXd b_l_new = (alpha.row(l).array() * mu.row(l).array()).matrix();
            b_bar = b_bar - b_l + b_l_new;
        }

        // Expected residual sum of squares, matching susieR's get_ER2():
        //   ERSS = ||y - X b_bar||^2 - sum_l ||X b_l||^2 + sum_l sum_j d_j alpha_lj mu2_lj
        // The exact ||X b_l||^2 (a full quadratic form) is required so that LD
        // cross-terms are retained; approximating it by sum_j d_j (alpha_lj mu_lj)^2
        // inflates sigma2 for correlated X and over-shrinks the single effects.
        VectorXd Xb = X * b_bar;
        double erss = (y - Xb).squaredNorm();
        for (int l = 0; l < L; ++l) {
            VectorXd b_l = (alpha.row(l).array() * mu.row(l).array()).matrix();
            erss -= (X * b_l).squaredNorm();
            ArrayXd a  = alpha.row(l).array();
            ArrayXd m2 = mu2.row(l).array();
            erss += (xtx * (a * m2)).sum();
        }
        if (opt.estimate_residual_variance && std::isfinite(erss) && erss > 0.0) sigma2 = erss / n;

        double elbo = -0.5 * n * std::log(2.0 * M_PI * sigma2) - erss / (2.0 * sigma2);
        for (int l = 0; l < L; ++l) elbo += lbf(l);
        fit.elbo.push_back(elbo);

        fit.niter = iter + 1;
        bool done = false;
        if (opt.convergence_method == SusieOptions::ELBO) {
            if (iter > 0 && std::abs(elbo - obj_prev) < opt.tol) done = true;
            obj_prev = elbo;
        } else {
            double da = (alpha - alpha_prev).cwiseAbs().maxCoeff();
            if (iter > 0 && da < opt.tol) done = true;
            alpha_prev = alpha;
        }
        if (done) { fit.converged = true; break; }
    }

    VectorXd pip(p);
    for (int j = 0; j < p; ++j) {
        double prod = 1.0;
        for (int l = 0; l < L; ++l) prod *= (1.0 - alpha(l, j));
        pip(j) = 1.0 - prod;
    }

    VectorXd Xb = X * b_bar;
    fit.alpha = alpha; fit.mu = mu; fit.mu2 = mu2;
    fit.lbf_variable = lbf_variable; fit.lbf = lbf;
    fit.pip = pip; fit.Xr = Xb; fit.fitted = Xb;
    fit.V = V; fit.sigma2 = sigma2; fit.intercept = 0.0;
    return fit;
}

// ---------------------------------------------------------------------------
// SuSiE-inf (unmappable infinitesimal effects), matching susieR's
// unmappable_effects = "inf" with estimate_residual_method = "MoM" and
// PIP-based convergence.
//
// Model: y = sum_l X b_l + X theta + e, theta_j ~ N(0, tau2), e_i ~ N(0, sigma2).
// The infinitesimal effect makes the residual covariance tau2 XX' + sigma2 I.
// All per-effect single-effect regressions are performed in the eigenspace of
// the (standardized) design so that Omega = (tau2 XX' + sigma2 I)^{-1} enters
// through Omega-weighted sufficient statistics X'Omega y and diag(X'Omega X).
// ---------------------------------------------------------------------------

// Negative SER log-likelihood as a function of the prior variance V, on the
// Omega-whitened scale (predictor weights pw = diag(X'Omega X), residual r =
// X'Omega (y - X b_{-l})). Mirrors susieR's neg_loglik on the inf path with the
// shat2 inflation factor equal to 1 (no finite-reference-R correction).
static double susie_inf_negll(double V, const ArrayXd& pw, const ArrayXd& res,
                              const ArrayXd& logpi) {
    ArrayXd denom = 1.0 + V * pw;
    ArrayXd lbf = -0.5 * denom.log() + 0.5 * V * res.square() / denom;
    ArrayXd w = lbf + logpi;
    double m = w.maxCoeff();
    return -(m + std::log((w - m).exp().sum()));
}

// Optimize V over [0,1] (golden section) then keep whichever of {optimum, V_init}
// has the higher likelihood, and finally snap to 0 when the null (V=0, loglik 0)
// is at least as good. Mirrors susieR's optimize_scalar_prior_variance for inf.
static double susie_inf_optimize_V(const ArrayXd& pw, const ArrayXd& res,
                                   const ArrayXd& logpi, double V_init) {
    const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
    double a = 0.0, b = 1.0;
    double c = b - gr * (b - a), d = a + gr * (b - a);
    double fc = susie_inf_negll(c, pw, res, logpi);
    double fd = susie_inf_negll(d, pw, res, logpi);
    for (int i = 0; i < 100 && (b - a) > 1e-6; ++i) {
        if (fc < fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = susie_inf_negll(c, pw, res, logpi); }
        else         { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = susie_inf_negll(d, pw, res, logpi); }
    }
    double v_opt = (a + b) / 2.0;
    double f_opt = susie_inf_negll(v_opt, pw, res, logpi);
    double f_init = susie_inf_negll(V_init, pw, res, logpi);
    double v_best = (f_init < f_opt) ? V_init : v_opt;
    double f_best = (f_init < f_opt) ? f_init : f_opt;
    // null (V=0) has loglik 0; keep it if it is at least as good
    if (f_best >= 0.0) return 0.0;
    return v_best;
}

SusieFit fit_susie_inf(const MatrixXd& X, const VectorXd& y, const SusieOptions& opt) {
    const int n = (int)X.rows();
    const int p = (int)X.cols();
    const int L = std::max(1, std::min(opt.L, p));
    const double var_y = y.squaredNorm() / std::max(1, n - 1); // y is pre-centered

    // --- thin eigendecomposition of standardized X via the n x n Gram matrix ---
    // (n < p in fine-mapping, so this is far cheaper than a full SVD). Components
    // with ~0 singular value contribute ~0 everywhere and are dropped, which
    // matches susieR's svd()-based path to numerical precision.
    MatrixXd G = X * X.transpose();                 // n x n
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(G);
    const VectorXd& evals = es.eigenvalues();       // ascending
    const MatrixXd& U_all = es.eigenvectors();
    const double ev_tol = std::max(evals.maxCoeff(), 0.0) * 1e-8;
    std::vector<int> keep;
    for (int k = 0; k < n; ++k) if (evals(k) > ev_tol) keep.push_back(k);
    const int r = (int)keep.size();

    VectorXd eigval(r);        // d_k^2
    MatrixXd Vmat(p, r);       // variant-space eigenvectors (like susieR's svd$v)
    VectorXd VtXty(r);         // d_k * (U_k' y) = V' X' y
    for (int idx = 0; idx < r; ++idx) {
        const int k = keep[idx];
        const double ev = evals(k);
        const double dk = std::sqrt(ev);
        VectorXd uk = U_all.col(k);
        eigval(idx)  = ev;
        Vmat.col(idx) = (X.transpose() * uk) / dk;
        VtXty(idx)   = dk * uk.dot(y);
    }
    MatrixXd Vsq = Vmat.array().square();           // p x r
    VectorXd Xty = X.transpose() * y;               // p
    const double yty = y.squaredNorm();
    ArrayXd logpi = ArrayXd::Constant(p, -std::log((double)p));

    // --- initialization (matches susieR init: uniform alpha, V = spv*var_y) ---
    MatrixXd alpha = MatrixXd::Constant(L, p, 1.0 / p);
    MatrixXd mu    = MatrixXd::Zero(L, p);
    MatrixXd mu2   = MatrixXd::Zero(L, p);
    MatrixXd lbf_variable = MatrixXd::Zero(L, p);
    VectorXd lbf   = VectorXd::Zero(L);
    VectorXd V     = VectorXd::Constant(L, opt.scaled_prior_variance * var_y);
    double sigma2  = var_y;
    double tau2    = 0.0;
    VectorXd theta = VectorXd::Zero(p);

    // Omega-derived caches, refreshed whenever (tau2, sigma2) change.
    ArrayXd omega_var = tau2 * eigval.array() + sigma2;                  // r
    VectorXd pw = Vsq * (eigval.array() / omega_var).matrix();          // p = diag(X'Omega X)
    VectorXd XtOmegay = Vmat * (VtXty.array() / omega_var).matrix();    // p

    VectorXd pip_prev = VectorXd::Zero(p);
    MatrixXd alpha_prev = alpha;

    SusieFit fit; fit.converged = false; fit.niter = 0;

    for (int iter = 0; iter < opt.max_iter; ++iter) {
        // ---- single-effect regressions (Omega-weighted) ----
        for (int l = 0; l < L; ++l) {
            VectorXd b_full  = (alpha.array() * mu.array()).colwise().sum().transpose();
            VectorXd b_l     = (alpha.row(l).array() * mu.row(l).array()).matrix().transpose();
            VectorXd b_minus = b_full - b_l;
            VectorXd Vtb     = Vmat.transpose() * b_minus;                          // r
            VectorXd XtOmegaXb = Vmat * (Vtb.array() * eigval.array() / omega_var).matrix();
            ArrayXd res = (XtOmegay - XtOmegaXb).array();                           // p (= X'Omega r)

            double Vl = susie_inf_optimize_V(pw.array(), res, logpi, V(l));
            V(l) = Vl;

            if (Vl <= 0.0) {
                alpha.row(l).setConstant(1.0 / p);
                mu.row(l).setZero();
                mu2.row(l).setZero();
                lbf_variable.row(l).setZero();
                lbf(l) = 0.0;
                continue;
            }
            ArrayXd denom = 1.0 + Vl * pw.array();
            ArrayXd lbfj  = -0.5 * denom.log() + 0.5 * Vl * res.square() / denom;
            ArrayXd w = lbfj + logpi;
            double m = w.maxCoeff();
            ArrayXd a = (w - m).exp();
            double s = a.sum();
            a /= s;
            ArrayXd post_var = Vl / denom;            // V*shat2/(V+shat2), shat2 = 1/pw
            ArrayXd post_mean = post_var * res;       // post_var/shat2 * betahat
            alpha.row(l) = a.transpose();
            mu.row(l)    = post_mean.transpose();
            mu2.row(l)   = (post_var + post_mean.square()).transpose();
            lbf_variable.row(l) = lbfj.transpose();
            lbf(l) = m + std::log(s);
        }
        fit.niter = iter + 1;

        // ---- PIP-based convergence (susieR forces this for inf) ----
        VectorXd pip(p);
        for (int j = 0; j < p; ++j) {
            double prod = 1.0;
            for (int l = 0; l < L; ++l) prod *= (1.0 - alpha(l, j));
            pip(j) = 1.0 - prod;
        }
        if (iter > 0) {
            double da = (alpha - alpha_prev).cwiseAbs().maxCoeff();
            double dp = (pip - pip_prev).cwiseAbs().maxCoeff();
            if (std::max(da, dp) < opt.tol) { fit.converged = true; break; }
        }
        alpha_prev = alpha;
        pip_prev   = pip;

        // ---- update variance components (method of moments) + theta BLUP ----
        VectorXd b   = (alpha.array() * mu.array()).colwise().sum().transpose();
        VectorXd Vtb_all = Vmat.transpose() * b;                                    // r
        // theta uses the current (pre-update) Omega caches and tau2, as in susieR.
        VectorXd XtOmegaXb_all = Vmat * (Vtb_all.array() * eigval.array() / omega_var).matrix();
        theta = tau2 * (XtOmegay - XtOmegaXb_all);

        // diag(V' M V), M = bbar bbar' - sum_l bl bl' + diag(sum_l alpha_l (mu_l^2 + 1/omega_l))
        ArrayXd diagVtMV = Vtb_all.array().square();
        ArrayXd tmpD = ArrayXd::Zero(p);
        for (int l = 0; l < L; ++l) {
            VectorXd bl   = (alpha.row(l).array() * mu.row(l).array()).matrix().transpose();
            VectorXd Vtbl = Vmat.transpose() * bl;
            diagVtMV -= Vtbl.array().square();
            ArrayXd omega_l = pw.array() + 1.0 / V(l);   // V(l)==0 -> +inf -> 1/omega_l == 0
            tmpD += alpha.row(l).transpose().array() *
                    (mu.row(l).transpose().array().square() + 1.0 / omega_l);
        }
        diagVtMV += (Vsq.transpose() * tmpD.matrix()).array();

        const double sumEv  = eigval.sum();
        const double sumEv2 = eigval.array().square().sum();
        Eigen::Matrix2d A;
        A << (double)n, sumEv, sumEv, sumEv2;
        double x1 = yty - 2.0 * b.dot(Xty) + (eigval.array() * diagVtMV).sum();
        double x2 = Xty.squaredNorm()
                    - 2.0 * (Vtb_all.array() * VtXty.array() * eigval.array()).sum()
                    + (eigval.array().square() * diagVtMV).sum();
        Eigen::Vector2d sol = A.colPivHouseholderQr().solve(Eigen::Vector2d(x1, x2));
        if (sol(0) > 0.0 && sol(1) > 0.0) { sigma2 = sol(0); tau2 = sol(1); }
        else                              { sigma2 = x1 / n; tau2 = 0.0; }
        if (sigma2 <= 0.0) sigma2 = var_y * 1e-8;

        // refresh Omega caches with the new (tau2, sigma2)
        omega_var = tau2 * eigval.array() + sigma2;
        pw        = Vsq * (eigval.array() / omega_var).matrix();
        XtOmegay  = Vmat * (VtXty.array() / omega_var).matrix();
    }

    VectorXd pip(p);
    for (int j = 0; j < p; ++j) {
        double prod = 1.0;
        for (int l = 0; l < L; ++l) prod *= (1.0 - alpha(l, j));
        pip(j) = 1.0 - prod;
    }

    VectorXd b_final = (alpha.array() * mu.array()).colwise().sum().transpose();
    VectorXd Xr = X * (b_final + theta);
    fit.alpha = alpha; fit.mu = mu; fit.mu2 = mu2;
    fit.lbf_variable = lbf_variable; fit.lbf = lbf;
    fit.pip = pip; fit.Xr = Xr; fit.fitted = Xr;
    fit.V = V; fit.sigma2 = sigma2; fit.intercept = 0.0;
    fit.tau2 = tau2; fit.theta = theta;
    return fit;
}

// ---------------------------------------------------------------------------
// SuSiE-ash (simplified port of SuSiE 2.0's unmappable_effects="ash").
//
// Model: y = sum_l X b_l + X theta + e,  e ~ N(0, sigma2 I),
//        theta_j ~ sum_k pi_k * N(0, sigma2 * sa2_k)      (sa2_0 = 0 = null)
//
// Per the SuSiE 2.0 paper (McCreight et al. 2025, Algorithm 1) the SER updates
// are performed under the Omega-weighted marginal likelihood, identical to the
// SuSiE-inf machinery, and tau2 = sigma2 * sum_k pi_k * sa2_k so the mixture
// prior implicitly defines Omega = (tau2 * X X' + sigma2 * I)^{-1}. Concretely,
// each outer iteration:
//   1. SER (Omega-weighted using cached omega_var / X'Omega y / diag(X'Omega X))
//   2. Provisional (sigma2, tau2) via MoM (2x2 solve, same as SuSiE-inf)
//   3. sa2 grid built from tau2/sigma2 so it spans plausible polygenic scales
//   4. Mr.ASH inner loop: coordinate-ascent NM updates for theta_j given the
//      residual r_theta_j = y - X*b_bar - X*theta + x_j*theta_j (primal scale)
//   5. pi = mean of gamma columns; tau2 = sigma2 * sum(pi * sa2)
//   6. Refresh Omega cache
// PIP-based convergence, matching susieR's ash path.
//
// Simplifications relative to susieR 0.16.5: no p x p LD matrix, no 3-tier
// masking layer, no Beta-Binomial slot-activity. These add ~440 R lines of
// stateful state that the paper says are for "elevated FDR under polygenic
// backgrounds" but are not required for the core algorithm.
// ---------------------------------------------------------------------------

// One sweep of Mr.ASH coordinate ascent. Updates theta, theta2 (E[theta_j^2]),
// and gamma (p x K responsibilities); r is updated incrementally so that at
// every point r = y - X*b_bar - X*theta. Returns max |delta theta_j|.
static double mr_ash_sweep(const MatrixXd& X, VectorXd& r,
                           VectorXd& theta, VectorXd& theta2, MatrixXd& gamma,
                           const ArrayXd& xtx, double sigma2,
                           const ArrayXd& sa2, const ArrayXd& log_pi) {
    const int p = (int)X.cols();
    const int K = (int)sa2.size();
    double max_delta = 0.0;
    for (int j = 0; j < p; ++j) {
        double xtx_j = xtx(j);
        if (!std::isfinite(xtx_j) || xtx_j <= 0.0) continue;

        double theta_old = theta(j);
        // partial residual: r_theta_j = r + x_j * theta_j_old
        r += theta_old * X.col(j);

        double xtr = X.col(j).dot(r);
        double bhat = xtr / xtx_j;

        // per-component NM stats on prior N(0, sigma2 * sa2_k):
        //   pvar_k  = sigma2 * sa2_k / (1 + sa2_k * xtx_j)
        //   pmean_k = (pvar_k / shat2) * bhat   with shat2 = sigma2 / xtx_j
        //   lbf_k   = -0.5*log(1 + sa2_k*xtx_j) + 0.5*bhat^2*sa2_k*xtx_j /
        //             (sigma2 * (1 + sa2_k*xtx_j))
        double m_max = -std::numeric_limits<double>::infinity();
        std::vector<double> logw(K), pmean(K), pvar(K);
        for (int k = 0; k < K; ++k) {
            double sk = sa2(k);
            if (sk <= 0.0) {
                pvar[k] = 0.0; pmean[k] = 0.0;
                logw[k] = log_pi(k);
            } else {
                double denom = 1.0 + sk * xtx_j;
                pvar[k]  = sigma2 * sk / denom;
                pmean[k] = sk * xtx_j / denom * bhat;    // == (pvar_k / shat2) * bhat
                double lbf = -0.5 * std::log(denom)
                           + 0.5 * bhat * bhat * sk * xtx_j / (sigma2 * denom);
                logw[k] = log_pi(k) + lbf;
            }
            if (logw[k] > m_max) m_max = logw[k];
        }
        double s = 0.0;
        std::vector<double> w(K);
        for (int k = 0; k < K; ++k) { w[k] = std::exp(logw[k] - m_max); s += w[k]; }
        double inv_s = 1.0 / s;

        double theta_new = 0.0, theta2_new = 0.0;
        for (int k = 0; k < K; ++k) {
            double gk = w[k] * inv_s;
            gamma(j, k) = gk;
            theta_new  += gk * pmean[k];
            theta2_new += gk * (pvar[k] + pmean[k] * pmean[k]);
        }
        theta(j)  = theta_new;
        theta2(j) = theta2_new;
        r -= theta_new * X.col(j);

        double d = std::abs(theta_new - theta_old);
        if (d > max_delta) max_delta = d;
    }
    return max_delta;
}

// Solve the SuSiE-inf 2x2 MoM system for (sigma2, tau2) and return the
// non-negative pair (or (x1/n, 0) when the system is not proper).
static void mom_variance_components(
    int n, double yty, const VectorXd& b, const VectorXd& Xty,
    const VectorXd& eigval, const VectorXd& Vtb_all, const VectorXd& VtXty,
    const ArrayXd& diagVtMV, double var_y_floor,
    double& sigma2, double& tau2) {
    const double sumEv  = eigval.sum();
    const double sumEv2 = eigval.array().square().sum();
    Eigen::Matrix2d A;
    A << (double)n, sumEv, sumEv, sumEv2;
    double x1 = yty - 2.0 * b.dot(Xty) + (eigval.array() * diagVtMV).sum();
    double x2 = Xty.squaredNorm()
                - 2.0 * (Vtb_all.array() * VtXty.array() * eigval.array()).sum()
                + (eigval.array().square() * diagVtMV).sum();
    Eigen::Vector2d sol = A.colPivHouseholderQr().solve(Eigen::Vector2d(x1, x2));
    if (sol(0) > 0.0 && sol(1) > 0.0) { sigma2 = sol(0); tau2 = sol(1); }
    else                              { sigma2 = std::max(x1 / n, var_y_floor); tau2 = 0.0; }
}

SusieFit fit_susie_ash(const MatrixXd& X, const VectorXd& y, const SusieOptions& opt) {
    const int n = (int)X.rows();
    const int p = (int)X.cols();
    const int L = std::max(1, std::min(opt.L, p));
    const bool fix_pi = !opt.ash_fix_pi.empty();
    const int K = fix_pi ? (int)opt.ash_fix_pi.size() : std::max(2, opt.ash_K);
    const double var_y = y.squaredNorm() / std::max(1, n - 1);

    // --- thin eigendecomposition of standardized X (same as fit_susie_inf) ---
    MatrixXd G = X * X.transpose();
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(G);
    const VectorXd& evals = es.eigenvalues();
    const MatrixXd& U_all = es.eigenvectors();
    const double ev_tol = std::max(evals.maxCoeff(), 0.0) * 1e-8;
    std::vector<int> keep;
    for (int k = 0; k < n; ++k) if (evals(k) > ev_tol) keep.push_back(k);
    const int r_rank = (int)keep.size();

    VectorXd eigval(r_rank);
    MatrixXd Vmat(p, r_rank);
    VectorXd VtXty(r_rank);
    for (int idx = 0; idx < r_rank; ++idx) {
        const int k = keep[idx];
        const double ev = evals(k);
        const double dk = std::sqrt(ev);
        VectorXd uk = U_all.col(k);
        eigval(idx)  = ev;
        Vmat.col(idx) = (X.transpose() * uk) / dk;
        VtXty(idx)   = dk * uk.dot(y);
    }
    MatrixXd Vsq = Vmat.array().square();
    VectorXd Xty = X.transpose() * y;
    const double yty = y.squaredNorm();
    ArrayXd logpi_ser = ArrayXd::Constant(p, -std::log((double)p));
    ArrayXd xtx = X.array().square().colwise().sum();
    for (int j = 0; j < p; ++j) if (xtx(j) <= 0.0) xtx(j) = std::numeric_limits<double>::infinity();

    // --- SER state ---
    MatrixXd alpha = MatrixXd::Constant(L, p, 1.0 / p);
    MatrixXd mu    = MatrixXd::Zero(L, p);
    MatrixXd mu2   = MatrixXd::Zero(L, p);
    MatrixXd lbf_variable = MatrixXd::Zero(L, p);
    VectorXd lbf   = VectorXd::Zero(L);
    VectorXd V     = VectorXd::Constant(L, opt.scaled_prior_variance * var_y);
    double sigma2  = var_y;
    double tau2    = 0.0;

    // --- Mr.ASH state ---
    // Grid: sa2_0 = 0 (null); non-null components log-spaced. The scale of the
    // grid is rebuilt each outer iteration from the MoM tau2/sigma2 so it stays
    // matched to the polygenic magnitude the data actually implies.
    ArrayXd sa2 = ArrayXd::Zero(K);
    const bool fix_sa2 = !opt.ash_fix_sa2.empty();
    if (fix_sa2) {
        // If not co-set with fix_pi, K falls back to ash_K; caller must ensure
        // ash_fix_sa2.size() == K. Silently clip to min length otherwise.
        int nsa2 = (int)opt.ash_fix_sa2.size();
        for (int k = 0; k < K && k < nsa2; ++k) sa2(k) = opt.ash_fix_sa2[k];
    } else {
        double sa2_max = 1.0;                       // initial guess ~ N(0,sigma2)
        double sa2_min = sa2_max * 1e-3;
        double lo = std::log(sa2_min), hi = std::log(sa2_max);
        sa2(0) = 0.0;
        for (int k = 1; k < K; ++k)
            sa2(k) = std::exp(lo + (hi - lo) * (double)(k - 1) / (double)std::max(1, K - 2));
    }
    ArrayXd mix_pi = ArrayXd::Constant(K, 1.0 / (double)K);
    if (fix_pi) {
        for (int k = 0; k < K; ++k) mix_pi(k) = opt.ash_fix_pi[k];
        double s = mix_pi.sum();
        if (s > 0.0) mix_pi /= s;
    } else {
        // Init pi with the null component dominant so ash does not absorb signal
        // before SER has had a chance to fine-map it.
        mix_pi(0) = 0.9; for (int k = 1; k < K; ++k) mix_pi(k) = 0.1 / (double)(K - 1);
    }
    ArrayXd log_pi_arr = mix_pi.max(1e-300).log();

    VectorXd theta  = VectorXd::Zero(p);
    VectorXd theta2 = VectorXd::Zero(p);
    MatrixXd gamma  = MatrixXd::Zero(p, K); gamma.col(0).setOnes();

    // Omega caches (theta enters through Omega, not by subtracting X*theta from y).
    ArrayXd omega_var = tau2 * eigval.array() + sigma2;
    VectorXd pw = Vsq * (eigval.array() / omega_var).matrix();
    VectorXd XtOmegay = Vmat * (VtXty.array() / omega_var).matrix();

    // maintained primal residual for the Mr.ASH inner loop: r = y - X*b_bar - X*theta
    VectorXd b_bar = VectorXd::Zero(p);
    VectorXd r_primal = y;

    MatrixXd alpha_prev = alpha;
    VectorXd pip_prev   = VectorXd::Zero(p);
    SusieFit fit; fit.converged = false; fit.niter = 0;

    for (int iter = 0; iter < opt.max_iter; ++iter) {
        // ---- Step 1: Omega-weighted single-effect regressions ----
        for (int l = 0; l < L; ++l) {
            VectorXd b_full  = (alpha.array() * mu.array()).colwise().sum().transpose();
            VectorXd b_l     = (alpha.row(l).array() * mu.row(l).array()).matrix().transpose();
            VectorXd b_minus = b_full - b_l;
            VectorXd Vtb     = Vmat.transpose() * b_minus;
            VectorXd XtOmegaXb = Vmat * (Vtb.array() * eigval.array() / omega_var).matrix();
            ArrayXd res = (XtOmegay - XtOmegaXb).array();

            double Vl = susie_inf_optimize_V(pw.array(), res, logpi_ser, V(l));
            V(l) = Vl;
            if (Vl <= 0.0) {
                alpha.row(l).setConstant(1.0 / p);
                mu.row(l).setZero();
                mu2.row(l).setZero();
                lbf_variable.row(l).setZero();
                lbf(l) = 0.0;
                continue;
            }
            ArrayXd denom = 1.0 + Vl * pw.array();
            ArrayXd lbfj  = -0.5 * denom.log() + 0.5 * Vl * res.square() / denom;
            ArrayXd w = lbfj + logpi_ser;
            double m = w.maxCoeff();
            ArrayXd a = (w - m).exp();
            double s = a.sum();
            a /= s;
            ArrayXd post_var = Vl / denom;
            ArrayXd post_mean = post_var * res;
            alpha.row(l) = a.transpose();
            mu.row(l)    = post_mean.transpose();
            mu2.row(l)   = (post_var + post_mean.square()).transpose();
            lbf_variable.row(l) = lbfj.transpose();
            lbf(l) = m + std::log(s);
        }
        fit.niter = iter + 1;

        // ---- PIP-based convergence ----
        VectorXd pip(p);
        for (int j = 0; j < p; ++j) {
            double prod = 1.0;
            for (int l = 0; l < L; ++l) prod *= (1.0 - alpha(l, j));
            pip(j) = 1.0 - prod;
        }
        if (iter > 0) {
            double da = (alpha - alpha_prev).cwiseAbs().maxCoeff();
            double dp = (pip - pip_prev).cwiseAbs().maxCoeff();
            if (std::max(da, dp) < opt.tol) { fit.converged = true; break; }
        }
        alpha_prev = alpha;
        pip_prev   = pip;

        // ---- Step 2: provisional (sigma2, tau2) via MoM, same 2x2 as SuSiE-inf ----
        VectorXd b = (alpha.array() * mu.array()).colwise().sum().transpose();
        VectorXd Vtb_all = Vmat.transpose() * b;
        ArrayXd diagVtMV = Vtb_all.array().square();
        ArrayXd tmpD = ArrayXd::Zero(p);
        for (int l = 0; l < L; ++l) {
            VectorXd bl   = (alpha.row(l).array() * mu.row(l).array()).matrix().transpose();
            VectorXd Vtbl = Vmat.transpose() * bl;
            diagVtMV -= Vtbl.array().square();
            ArrayXd omega_l = pw.array() + 1.0 / V(l);
            tmpD += alpha.row(l).transpose().array() *
                    (mu.row(l).transpose().array().square() + 1.0 / omega_l);
        }
        diagVtMV += (Vsq.transpose() * tmpD.matrix()).array();

        double sigma2_prov = sigma2, tau2_prov = tau2;
        mom_variance_components(n, yty, b, Xty, eigval, Vtb_all, VtXty, diagVtMV,
                                 var_y * 1e-8, sigma2_prov, tau2_prov);
        if (opt.estimate_residual_variance && sigma2_prov > 0.0) sigma2 = sigma2_prov;

        // ---- Step 3: build a data-driven sa2 grid based on tau2_prov/sigma2 ----
        // Interpret tau2 as an average variance sigma2*E[sa2]; set sa2_max so a
        // few non-null components cover several times this magnitude, and
        // sa2_min three orders of magnitude below, log-spaced.
        // Skip grid rebuild when pi or sa2 is fixed (reduction tests want a pinned grid).
        if (!fix_pi && !fix_sa2) {
            double ratio = (sigma2 > 0.0) ? std::max(tau2_prov / sigma2, 1e-8) : 1e-4;
            double sa2_max = std::max(4.0 * ratio, 1e-3);
            double sa2_min = std::max(sa2_max * 1e-3, 1e-8);
            double lo = std::log(sa2_min), hi = std::log(sa2_max);
            sa2(0) = 0.0;
            for (int k = 1; k < K; ++k)
                sa2(k) = std::exp(lo + (hi - lo) * (double)(k - 1) / (double)std::max(1, K - 2));
        }

        // ---- Step 4: Mr.ASH inner loop for theta ----
        // Sync r_primal to the fresh b_bar before sweeping.
        VectorXd b_bar_new = b;
        r_primal += X * (b_bar - b_bar_new);          // undo old b_bar, apply new
        b_bar = b_bar_new;
        // r_primal now = y - X*b_bar - X*theta_old.

        double inner_tol = std::max(opt.ash_inner_tol, opt.tol);
        for (int inner = 0; inner < opt.ash_inner_max; ++inner) {
            double d = mr_ash_sweep(X, r_primal, theta, theta2, gamma,
                                     xtx, sigma2, sa2, log_pi_arr);
            if (d < inner_tol) break;
        }

        // ---- Step 5: EM update of pi = mean of gamma columns, tau2 = sigma2 * sum(pi*sa2) ----
        if (!fix_pi) {
            Eigen::RowVectorXd cm = gamma.colwise().mean();
            mix_pi = cm.array().transpose();
            mix_pi = mix_pi.max(1e-8);
            mix_pi /= mix_pi.sum();
            log_pi_arr = mix_pi.log();
        }
        tau2 = sigma2 * (mix_pi * sa2).sum();

        // ---- Step 6: refresh Omega cache ----
        omega_var = tau2 * eigval.array() + sigma2;
        pw        = Vsq * (eigval.array() / omega_var).matrix();
        XtOmegay  = Vmat * (VtXty.array() / omega_var).matrix();
    }

    VectorXd pip(p);
    for (int j = 0; j < p; ++j) {
        double prod = 1.0;
        for (int l = 0; l < L; ++l) prod *= (1.0 - alpha(l, j));
        pip(j) = 1.0 - prod;
    }

    VectorXd b_final = (alpha.array() * mu.array()).colwise().sum().transpose();
    VectorXd Xr = X * (b_final + theta);
    fit.alpha = alpha; fit.mu = mu; fit.mu2 = mu2;
    fit.lbf_variable = lbf_variable; fit.lbf = lbf;
    fit.pip = pip; fit.Xr = Xr; fit.fitted = Xr;
    fit.V = V; fit.sigma2 = sigma2; fit.intercept = 0.0;
    fit.tau2 = tau2; fit.theta = theta;
    fit.ash_sa2 = sa2.matrix();
    fit.ash_pi  = mix_pi.matrix();
    return fit;
}

// ---------------------------------------------------------------------------
// Credible sets
// ---------------------------------------------------------------------------

// Minimum absolute correlation among a set of (standardized) columns of X_std.
// For large sets the estimate is based on a random subsample of up to 100 columns
// (mirrors susieR's get_purity()).
static double set_purity(const std::vector<int32_t>& vars, const MatrixXd& X_std) {
    const int k = (int)vars.size();
    if (k <= 1) return 1.0;

    std::vector<int32_t> use = vars;
    if (k > 100) {
        std::mt19937 rng(12345);
        std::shuffle(use.begin(), use.end(), rng);
        use.resize(100);
    }
    const int m = (int)use.size();
    const double denom = std::max(1.0, (double)(X_std.rows() - 1)); // ||col||^2 for unit-sd columns
    double min_abs = 1.0;
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            double num = X_std.col(use[i]).dot(X_std.col(use[j]));
            double denom_ij = std::sqrt(X_std.col(use[i]).squaredNorm() * X_std.col(use[j]).squaredNorm());
            double corr = (denom_ij > 0.0) ? (num / denom_ij) : 0.0;
            double a = std::abs(corr);
            if (a < min_abs) min_abs = a;
        }
    }
    (void)denom;
    return min_abs;
}

std::vector<CredibleSet> susie_get_cs(const SusieFit& fit, const MatrixXd& X_std,
                                      double coverage, double min_abs_corr) {
    std::vector<CredibleSet> out;
    const int L = (int)fit.alpha.rows();
    const int p = (int)fit.alpha.cols();
    if (L == 0 || p == 0) return out;

    std::vector<std::vector<int32_t> > kept_sets; // for de-duplication

    for (int l = 0; l < L; ++l) {
        // order variables by descending alpha within this effect
        std::vector<int32_t> order(p);
        std::iota(order.begin(), order.end(), 0);
        const double* arow = fit.alpha.row(l).data();
        std::sort(order.begin(), order.end(),
                  [&](int32_t a, int32_t b) { return fit.alpha(l, a) > fit.alpha(l, b); });

        // accumulate until cumulative alpha reaches the requested coverage
        double cum = 0.0;
        std::vector<int32_t> vars;
        for (int idx = 0; idx < p; ++idx) {
            int32_t j = order[idx];
            vars.push_back(j);
            cum += fit.alpha(l, j);
            if (cum >= coverage) break;
        }
        (void)arow;

        // an effect that explains nothing spreads alpha ~ uniformly; its set will be
        // large with low purity and is filtered out below.
        double purity = set_purity(vars, X_std);
        if (purity < min_abs_corr) continue;

        // de-duplicate: skip if an identical set was already reported
        std::vector<int32_t> sorted_vars = vars;
        std::sort(sorted_vars.begin(), sorted_vars.end());
        bool dup = false;
        for (const auto& ks : kept_sets) if (ks == sorted_vars) { dup = true; break; }
        if (dup) continue;
        kept_sets.push_back(sorted_vars);

        CredibleSet cs;
        cs.effect_index = l;
        cs.variables = vars;
        cs.pips.reserve(vars.size());
        for (int32_t j : vars) cs.pips.push_back(fit.pip(j));
        cs.coverage = cum;
        cs.requested_coverage = coverage;
        cs.purity = purity;
        cs.lbf = fit.lbf(l);
        cs.V = fit.V(l);
        out.push_back(std::move(cs));
    }
    return out;
}

// ---------------------------------------------------------------------------
// High-level wrapper
// ---------------------------------------------------------------------------

SusieResult simple_susie_without_missing(const VectorXd& pheno_vec,
                                         const MatrixXd& geno_mat,
                                         const SusieOptions& opt,
                                         double coverage,
                                         double min_abs_corr) {
    const int n = (int)geno_mat.rows();
    const int p = (int)geno_mat.cols();

    // center the phenotype (covariates already regressed out by the caller)
    VectorXd y = pheno_vec.array() - pheno_vec.mean();
    if (!y.allFinite())
        error("simple_susie_without_missing(): the phenotype vector contains non-finite values");

    // Work on a standardized copy of X so the caller's matrix is left untouched.
    //
    // Every single-effect regression shares one residual and one softmax across all p
    // variables, so a single bad column contaminates the entire fit rather than just
    // its own coefficient (unlike the marginal per-variant regressions, which stay
    // correct for every other variant). Two kinds of bad column are screened here:
    //   - non-finite entries (e.g. a NaN dosage that survived mean-imputation);
    //   - zero variance (monomorphic in the analyzed samples, an imputed variant with
    //     an identical dosage for everybody, or a variant fully explained by the
    //     covariates after the FWL residualization).
    // Both are zeroed out and reported; fit_susie() then gives them zero prior weight
    // so their PIP/alpha are 0 and the remaining variants are fit normally.
    MatrixXd X = geno_mat;
    int n_nonfinite_cols = 0, n_novar_cols = 0;
    std::vector<int> bad_cols;
    for (int j = 0; j < p; ++j) {
        bool finite = X.col(j).allFinite();
        double mean = finite ? X.col(j).mean() : 0.0;
        double ss   = finite ? (X.col(j).array() - mean).square().sum() : 0.0;
        double sd   = std::sqrt(ss / std::max(1, n - 1));
        if (!finite || !(ss > 0.0) || !(sd > 0.0)) {
            if (!finite) ++n_nonfinite_cols; else ++n_novar_cols;
            if ((int)bad_cols.size() < 10) bad_cols.push_back(j);
            X.col(j).setZero();
            continue;
        }
        if (opt.standardize) X.col(j) = (X.col(j).array() - mean) / sd;
        else                 X.col(j) = X.col(j).array() - mean; // at least center
    }
    if (n_nonfinite_cols + n_novar_cols > 0) {
        std::string idx;
        for (size_t i = 0; i < bad_cols.size(); ++i) {
            if (i) idx += ",";
            idx += std::to_string(bad_cols[i]);
        }
        if ((int)bad_cols.size() < n_nonfinite_cols + n_novar_cols) idx += ",...";
        notice("SuSiE: excluding %d of %d variants with no usable genotype variance "
               "(%d with non-finite values, %d monomorphic/collinear-with-covariates); "
               "0-based column indices: %s",
               n_nonfinite_cols + n_novar_cols, p, n_nonfinite_cols, n_novar_cols, idx.c_str());
    }

    SusieResult res;
    if (opt.unmappable_effects == SusieOptions::INF)
        res.fit = fit_susie_inf(X, y, opt);
    else if (opt.unmappable_effects == SusieOptions::ASH)
        res.fit = fit_susie_ash(X, y, opt);
    else
        res.fit = fit_susie(X, y, opt);
    res.cs = susie_get_cs(res.fit, X, coverage, min_abs_corr);
    return res;
}

} // namespace susie
