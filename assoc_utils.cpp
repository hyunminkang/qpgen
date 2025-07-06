#include "assoc_utils.h"
#include "qgenlib/qgen_error.h"

Eigen::MatrixXd bulk_adjust_for_covariates(const Eigen::MatrixXd& values, const Eigen::MatrixXd& covariates) {
    const long g = values.rows();
    const long n = values.cols();
    const long p = covariates.rows();

    // --- Input Validation ---
    // Ensure the number of observations (columns) is the same in both matrices.
    if (n != covariates.cols()) {
        throw std::runtime_error("Value and covariate matrices must have the same number of columns (observations).");
    }

    // --- Transpose Inputs for Regression ---
    // The standard regression setup is (observations x variables).
    // We transpose our g x n and p x n inputs to n x g and n x p.
    Eigen::MatrixXd V_t = values.transpose(); // V_t is n x g
    Eigen::MatrixXd C_t = covariates.transpose(); // C_t is n x p

    // --- Augment Covariate Matrix with Intercept ---
    // We create an augmented covariate matrix `C_aug` of size n x (p+1).
    Eigen::MatrixXd C_aug(n, p + 1);
    C_aug.col(0).setOnes(); // First column is the intercept (all ones).
    C_aug.rightCols(p) = C_t; // The remaining columns are the transposed covariates.

    // --- Solve the Least Squares Problem ---
    // We solve for the coefficient matrix `beta` ((p+1) x g) where `V_t = C_aug * beta`.
    // Eigen's `colPivHouseholderQr()` provides a robust way to solve this.
    Eigen::MatrixXd beta = C_aug.colPivHouseholderQr().solve(V_t);

    // --- Calculate Residuals ---
    // The predicted values are `C_aug * beta`.
    // The residuals are the original (transposed) values minus the predicted values.
    // The resulting residual matrix is in n x g format.
    Eigen::MatrixXd residuals_t = V_t - (C_aug * beta);

    // --- Transpose Result Back to g x n Format ---
    return residuals_t.transpose();
}

void center_rows(Eigen::MatrixXd &matrix) {
    // We iterate over each row of the matrix.
    // The .rowwise() method allows us to apply an operation to each row.
    // The .mean() method calculates the mean of the elements in a row.
    // This creates a column vector where each element is the mean of the corresponding row in the original matrix.
    Eigen::VectorXd row_means = matrix.rowwise().mean();

    notice("Centering rows of matrix with dimensions %d x %d.", matrix.rows(), matrix.cols());

    matrix.colwise() -= row_means;

    // // We then iterate over each row again to subtract the calculated mean.
    // for (int i = 0; i < matrix.rows(); ++i) {
    //     // .noalias() is a performance optimization in Eigen. It tells Eigen that
    //     // the result of the expression on the right-hand side does not alias
    //     // with the matrix on the left-hand side, allowing for a more efficient
    //     // subtraction operation without creating temporary objects.
    //     matrix.row(i).noalias() -= Eigen::VectorXd::Constant(matrix.cols(), row_means(i));
    // }
}

/**
 * @brief Computes the continued fraction part for the incomplete beta function.
 *
 * This is an implementation of Lentz's method as described in Numerical Recipes.
 * It's a key component for accurately calculating the incomplete beta function.
 *
 * @param a Parameter 'a' of the beta function.
 * @param b Parameter 'b' of the beta function.
 * @param x The upper limit of integration.
 * @return The result of the continued fraction evaluation.
 */
double betacf(double a, double b, double x) {
    const int MAXIT = 200;
    const double EPS = 1.0e-14;
    const double FPMIN = 1.0e-30;

    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= MAXIT; m++) {
        int m2 = 2 * m;
        // Even step
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        // Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) < EPS) {
            return h;
        }
    }
    // This warning can be enabled for debugging convergence issues.
    // std::cerr << "Warning: betacf did not converge for a=" << a << ", b=" << b << ", x=" << x << std::endl;
    return h;
}

/**
 * @brief Computes the natural logarithm of the regularized incomplete beta function I_x(a,b).
 *
 * This function computes log(I_x(a,b)) using a continued fraction evaluation (betacf).
 * It uses a symmetry property to maintain accuracy when x is large.
 * The formula used is:
 * log(I_x(a,b)) = a*log(x) + b*log(1-x) - log(a) - log(B(a,b)) + log(betacf(a,b,x))
 * where B(a,b) is the beta function.
 *
 * @param x The upper limit of integration (must be in [0, 1]).
 * @param a Parameter 'a' of the beta function (must be > 0).
 * @param b Parameter 'b' of the beta function (must be > 0).
 * @return The natural logarithm of I_x(a,b).
 */
double betaincln(double x, double a, double b) {
    if (x < 0.0 || x > 1.0) return std::numeric_limits<double>::quiet_NaN();
    if (x == 0.0) return -std::numeric_limits<double>::infinity();
    if (x == 1.0) return 0.0;

    // Use symmetry property for better accuracy when x is large
    if (x > (a + 1.0) / (a + b + 2.0)) {
        // Compute log(1 - I_{1-x}(b,a))
        // First compute log(I_{1-x}(b,a))
        double log_beta_func_sym = std::lgamma(b) + std::lgamma(a) - std::lgamma(b + a);
        double cf_sym = betacf(b, a, 1.0 - x);
        double log_power_term_sym = b * std::log(1.0 - x) + a * std::log(x);
        double log_i_sym = log_power_term_sym - std::log(b) - log_beta_func_sym + std::log(cf_sym);

        // log(1-y) is computed as log1p(-y) for precision when y is small.
        return std::log1p(-std::exp(log_i_sym));
    } else {
        // Direct computation for smaller x
        double log_beta_func = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
        double cf = betacf(a, b, x);
        double log_power_term = a * std::log(x) + b * std::log(1.0 - x);
        return log_power_term - std::log(a) - log_beta_func + std::log(cf);
    }
}

/**
 * @brief Converts a t-statistic to a 2-sided -log10 p-value.
 *
 * This function calculates the p-value for a given t-statistic and degrees
 * of freedom. It is designed to be numerically stable even for extremely
 * large t-statistic values by working in log-space to avoid underflow.
 * The p-value is derived from the student's t-distribution cumulative
 * distribution function, which is related to the regularized incomplete
 * beta function.
 *
 * @param tstat The t-statistic value.
 * @param df The degrees of freedom.
 * @return The 2-sided -log10 p-value. Returns NaN if df <= 0.
 */
double tstat2log10pval(double tstat, double df) {
    if (df <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double abs_t = std::abs(tstat);

    // The 2-sided p-value for the t-distribution is given by I_x(df/2, 1/2),
    // where x = df / (df + t^2).
    double x = df / (df + abs_t * abs_t);

    // To maintain numerical precision for very small p-values, we compute
    // the natural logarithm of the p-value directly.
    double log_p_value = betaincln(x, df / 2.0, 0.5);

    // Convert the natural log p-value to -log10 p-value.
    // -log10(p) = -ln(p) / ln(10)
    return -(log_p_value / std::log(10.0));
}
/**
 * @brief Performs simple linear regression for each column of a matrix against a vector.
 * * Assumes that the input vector y and each column of the matrix X are centered
 * (have a mean of zero).
 *
 * @param y An Eigen::VectorXd of size n (dependent variable).
 * @param X An Eigen::MatrixXd of size n x p (independent variables).
 * @return A std::vector of RegressionResult structs, one for each column of X.
 */
bool simple_linear_regression_without_missing(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, std::vector<slr_sumstat_t>& results) {
    int n = y.size();
    int p = X.cols();

    // for(int32_t i = 0; i < 5; ++i) {
    //     notice("y[%d] = %.5g, X[%d,0] = %.5g", i, y[i], i, X(i, 0));
    // }

    if (n != X.rows()) {
        error("Vector and Matrix dimensions do not match : (%d x %d) vs (%d x %d)", n, 1, X.rows(), X.cols());
    }
    
    // Degrees of freedom for the t-distribution
    // For simple linear regression (one predictor), it's n - 2
    int df = n - 2;
    if (df <= 0) {
        error("Not enough data points to perform regression (n-2 = %d must be > 0).", df);
    }

    results.clear();
    results.resize(p);

    // 1. Calculate squared norms for each column of X (x_i' * x_i)
    Eigen::VectorXd x_col_sq_norms = X.colwise().squaredNorm();

    // 2. Calculate X' * y
    Eigen::VectorXd xt_y = X.transpose() * y;

    // 3. Calculate betas for all columns
    // beta_i = (x_i' * y) / (x_i' * x_i)
    Eigen::ArrayXd betas = xt_y.array() / x_col_sq_norms.array();

    // 4. Calculate Sum of Squared Errors (SSE) for each regression
    // sse_i = y'y - (x_i'y)^2 / (x_i'x_i)
    double y_sq_norm = y.squaredNorm();
    Eigen::ArrayXd sse = y_sq_norm - (xt_y.array().square() / x_col_sq_norms.array());

    // 5. Calculate standard errors for all betas
    // se(beta_i) = sqrt( (sse_i / (n-2)) / (x_i'x_i) )
    Eigen::ArrayXd std_errors = (sse / (df * x_col_sq_norms.array())).sqrt();

    // 6. Calculate t-statistics for all betas
    Eigen::ArrayXd t_values = betas / std_errors;

    for (int i = 0; i < p; ++i) {
        slr_sumstat_t& ss = results[i];
        ss.beta = betas[i];
        ss.se = std_errors[i];
        ss.tstat = t_values[i];
        ss.log10p = tstat2log10pval(ss.tstat, df);
        ss.n_obs = n; // Store the number of observations
    }

    //error("stop");

    return true;
}

bool simple_linear_regression_with_missing(const Eigen::VectorXd& y,
                                            const Eigen::Vector<bool, Eigen::Dynamic>& y_mask,
                                            const Eigen::MatrixXd& X,
                                            const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& X_mask,
                                            std::vector<slr_sumstat_t>& results) {
    const int n_total = y.size();
    const int p = X.cols();

    // --- Input Validation ---
    if (n_total != X.rows() || n_total != y_mask.size() ||
        X.rows() != X_mask.rows() || X.cols() != X_mask.cols()) {
        error("Data and mask dimensions do not match.");
        return false;
    }

    results.clear();
    results.resize(p);

    // --- Main loop to iterate over each column of X ---
    for (int i = 0; i < p; ++i) {
        std::vector<double> y_complete;
        std::vector<double> x_complete;
        y_complete.reserve(n_total);
        x_complete.reserve(n_total);

        // 1. Pairwise deletion: Collect all pairs where mask is true
        for (int j = 0; j < n_total; ++j) {
            // THE CORE CHANGE IS HERE: Check the boolean mask instead of isnan()
            if (y_mask(j) && X_mask(j, i)) {
                y_complete.push_back(y(j));
                x_complete.push_back(X(j, i));
            }
        }

        const int n_complete = y_complete.size();
        slr_sumstat_t& ss = results[i];
        ss.n_obs = n_complete;

        // 2. Check for sufficient data to perform regression
        const int df = n_complete - 2;
        if (df <= 0) {
            ss.beta = std::numeric_limits<double>::quiet_NaN();
            ss.se = std::numeric_limits<double>::quiet_NaN();
            ss.tstat = std::numeric_limits<double>::quiet_NaN();
            ss.log10p = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // 3. Map the std::vectors to Eigen vectors for calculation (no copy)
        Eigen::Map<Eigen::VectorXd> y_vec(y_complete.data(), n_complete);
        Eigen::Map<Eigen::VectorXd> x_vec(x_complete.data(), n_complete);

        // 4. Perform regression calculations on the complete data
        const double x_sq_norm = x_vec.squaredNorm();
        const double y_sq_norm = y_vec.squaredNorm();
        const double xt_y = x_vec.dot(y_vec);

        if (x_sq_norm == 0) {
            ss.beta = std::numeric_limits<double>::quiet_NaN();
            ss.se = std::numeric_limits<double>::quiet_NaN();
            ss.tstat = std::numeric_limits<double>::quiet_NaN();
            ss.log10p = std::numeric_limits<double>::quiet_NaN();
            continue;
        }
        
        ss.beta = xt_y / x_sq_norm;
        
        const double sse = y_sq_norm - (xt_y * xt_y) / x_sq_norm;
        ss.se = std::sqrt((sse / df) / x_sq_norm);
        
        if (ss.se > 0) {
            ss.tstat = ss.beta / ss.se;
            ss.log10p = tstat2log10pval(ss.tstat, df);
        } else {
            ss.tstat = std::numeric_limits<double>::quiet_NaN();
            ss.log10p = std::numeric_limits<double>::quiet_NaN();
        }
    }

    return true;
}