#ifndef __ASSOC_UTILS_H
#define __ASSOC_UTILS_H

#include "Eigen/Dense"
#include "qgenlib/qgen_error.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"

Eigen::MatrixXd bulk_adjust_for_covariates(const Eigen::MatrixXd& values, const Eigen::MatrixXd& covariates);

Eigen::MatrixXd pheno_adj_cov_nxt_without_missing(const Eigen::MatrixXd& values, const Eigen::MatrixXd& covariates);

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

bool simple_rect_regression_without_missing(
    const Eigen::MatrixXd& Y,
    const Eigen::MatrixXd& X,
    std::vector<std::vector<slr_sumstat_t> >& results);

Eigen::VectorXd rint_with_missing(
    const Eigen::VectorXd& values,
    const Eigen::Vector<bool, Eigen::Dynamic>& mask);

Eigen::VectorXd rint_without_missing(const Eigen::VectorXd& values);

Eigen::MatrixXd rint_matrix_without_missing(const Eigen::MatrixXd& matrix);

void standardize_matrix_columns_inplace(Eigen::MatrixXd& matrix);

Eigen::VectorXd columnwise_dot(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2);

int32_t assoc_single_trait( 
    htsFile* wf, // output file handle
    const char* pheno_id, // phenotype ID
    const Eigen::VectorXd& phe_vec,      // phenotype vector
    const Eigen::VectorXd& phe_rint_vec, // rinted phenotype vector 
    const Eigen::Vector<bool, Eigen::Dynamic>& phe_mask_vec, // phenotype mask vector
    bool phe_has_missing, // if the phenotype has missing values
    bool skip_rint, // if the rank-based inverse normal transformation should be skipped
    const Eigen::MatrixXd& geno_mat,    // genotype matrix
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& geno_mask,   // genotype mask matrix
    bool geno_has_missing, // if the genotype has missing values
    const std::vector<cpra_t>& v_cpra,      // variant pairs
    const std::vector<int32_t>& ans,         // allele counts
    const std::vector<double>& acs,         // allele counts
    const std::vector<double>& infos       // infor values
);

int32_t assoc_single_trait( 
    htsFile* wf, // output file handle
    const char* pheno_id, // phenotype ID
    const Eigen::VectorXd& phe_vec,      // phenotype vector
    const Eigen::VectorXd& phe_rint_vec, // rinted phenotype vector 
    const Eigen::Vector<bool, Eigen::Dynamic>& phe_mask_vec, // phenotype mask vector
    bool phe_has_missing, // if the phenotype has missing values
    bool skip_rint, // if the rank-based inverse normal transformation should be skipped
    const Eigen::MatrixXd& geno_mat,    // genotype matrix
    const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& geno_mask,   // genotype mask matrix
    bool geno_has_missing, // if the genotype has missing values
    const std::vector<cpra_t>& v_cpra,      // variant pairs
    const std::vector<int32_t>& ans,         // allele counts
    const std::vector<double>& acs,         // allele counts
    const std::vector<int32_t>& gc0s,        // genotype counts for genotype 0
    const std::vector<int32_t>& gc1s,        // genotype counts for genotype 1
    const std::vector<int32_t>& gc2s,        // genotype counts for genotype 2
    const std::vector<double>& infos       // infor values
);

#endif // __ASSOC_UTILS_H