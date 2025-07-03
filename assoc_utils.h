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
};

typedef struct _slr_sumstat_t slr_sumstat_t;

double tstat2log10pval(double tstat, double df);
bool simple_linear_regression(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, std::vector<slr_sumstat_t>& results);

#endif // __ASSOC_UTILS_H