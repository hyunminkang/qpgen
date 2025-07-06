#ifndef __ASSOC_UTILS_H
#define __ASSOC_UTILS_H

#include "Eigen/Dense"

Eigen::MatrixXd bulk_adjust_for_covariates(const Eigen::MatrixXd& values, const Eigen::MatrixXd& covariates);

void center_rows(Eigen::MatrixXd &matrix);

// Structure to hold the results of the regression for one variable
struct _slr_sumstat_t {
    double beta;
    double se;
    double tstat;
    double log10p;
    int32_t n_obs;
};

typedef struct _slr_sumstat_t slr_sumstat_t;

double tstat2log10pval(double tstat, double df);
bool simple_linear_regression_without_missing(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, std::vector<slr_sumstat_t>& results);

bool simple_linear_regression_with_missing( const Eigen::VectorXd& y,
                                            const Eigen::Vector<bool, Eigen::Dynamic>& y_mask,
                                            const Eigen::MatrixXd& X,
                                            const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& X_mask,
                                            std::vector<slr_sumstat_t>& results);

#endif // __ASSOC_UTILS_H