#ifndef __ASSOC_UTILS_H
#define __ASSOC_UTILS_H

#include "Eigen/Dense"
#include "qgenlib/qgen_error.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"
#include "qpgen_utils.h"
#include "pheno.h"


class genotype_chunk {
public:
    Eigen::MatrixXd geno_mat; // genotype matrix (n_samples x n_variants)
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> geno_mask; // genotype missingness mask (n_samples x n_variants)
    std::vector<cpra_t> v_cpra; // vector of CPRA identifiers for the variants in the chunk
    std::vector<double> acs; // allele counts for the variants in the chunk
    std::vector<int32_t> ans; // allele numbers for the variants in the chunk
    std::vector<int32_t> gc0s, gc1s, gc2s;
    std::vector<double> infos;
    int32_t n_variants; // number of variants in the chunk
    int32_t n_skipped;
    int32_t n_samples; // number of samples in the chunk

    genotype_chunk() : n_variants(0), n_skipped(0), n_samples(0) {}
    ~genotype_chunk() {}
};

// class containing individual-level association analysis input data
class ind_assoc_input {
protected:
    bool load_pgenlist(const char* pgenlistf);
    bool load_pgen_files(const char* pgenf, const char* pivarf, const char* psamf);
    bool load_pheno_cov_matrices(const char* phef, const char* pheno_format, const char* covf, const char* cov_format);
    bool is_pheno_loaded() const { return !pheno_matrix.samp_ids.empty() && !pheno_matrix.pheno_ids.empty(); }
    //bool add_genotypes();

public:
    MultiPgenIdxReader mpr; // multi-pgen reader
    PhenoMatrix pheno_matrix; // phenotype matrix
    PhenoMatrix cov_matrix;   // covariate matrix
    std::vector<std::string> overlapping_sample_ids; // overlapping sample IDs among genotype, phenotype, and covariate files
    int32_t jump_thres_bp; // jump threshold in base pairs for the multi-pgen reader
    int32_t icol_pivar_idx; // column index for the variant ID in the pivar file
    std::vector<std::string> subset_sample_ids;
    std::vector<std::string> subset_pheno_ids;
    bool rint_before_adj; // Perform rank-based inverse normal transformation before covariate adjustment
    bool rint_after_adj;  // Perform rank-based inverse normal transformation after covariate adjustment
    bool geno_chunk_done; // true once the current region has been fully streamed
    double min_af; // minimum allele frequency for filtering variants
    double max_af; // maximum allele frequency for filtering variants
    double min_ac; // minimum allele count for filtering variants
    double max_ac; // maximum allele count for filtering variants
    genotype_chunk geno_chunk; // genotype chunk for the current region

    ind_assoc_input() {
        jump_thres_bp = 1000000;
        icol_pivar_idx = 8;
        rint_before_adj = false;
        rint_after_adj = false;
        geno_chunk_done = false;
    }
    ~ind_assoc_input() {}

    void set_jump_thres_bp(int32_t thres) { 
        jump_thres_bp = thres; 
        mpr.set_jump_thres_bp(thres);
    }
    void set_icol_pivar_idx(int32_t idx) {
        icol_pivar_idx = idx;
        mpr.set_icol_pivar_idx(idx);
    }
    void set_rint_before_adj(bool val) { rint_before_adj = val; }
    void set_rint_after_adj(bool val) { rint_after_adj = val; }

    void set_minmax_af(double _min_af, double _max_af) {
        min_af = _min_af;
        max_af = _max_af;
    }

    void set_minmax_ac(double _min_ac, double _max_ac) {
        min_ac = _min_ac;
        max_ac = _max_ac;
    }

    inline void set_subset_sample_file(const char* samplef) {
        if ( is_pheno_loaded() ) {
            error("Cannot set subset sample file after the phenotype matrix is loaded");
        }
        tsv_reader tr_sample(samplef);
        while(tr_sample.read_line()) {
            subset_sample_ids.push_back(tr_sample.str_field_at(0));
        }
    }

    inline void set_subset_sample_ids(const std::vector<std::string>& samp_ids) {
        if ( is_pheno_loaded() ) {
            error("Cannot set subset sample file after the phenotype matrix is loaded");
        }
        subset_sample_ids = samp_ids;
    }

    inline void set_subset_pheno_file(const char* phenof, int32_t icol_pheno_id = 1) {
        if ( is_pheno_loaded() ) {
            error("Cannot set subset phenotype file after the phenotype matrix is loaded");
        }
        tsv_reader tr_pheno(phenof);
        while(tr_pheno.read_line()) {
            subset_pheno_ids.push_back(tr_pheno.str_field_at(icol_pheno_id - 1));
        }
    }

    inline void set_subset_pheno_ids(const std::vector<std::string>& pheno_ids) {
        if ( is_pheno_loaded() ) {
            error("Cannot set subset phenotype file after the phenotype matrix is loaded");
        }
        subset_pheno_ids = pheno_ids;
    }

    bool process_pgenlist(const char* pgenlistf, const char* phef, const char* pheno_format, const char* covf, const char* cov_format);
    bool process_single_pgen(const char* pgenf, const char* pivarf, const char* psamf, const char* phef, const char* pheno_format, const char* covf, const char* cov_format);
    bool load_genotype_chunk(const char* chrom, int32_t beg, int32_t end, int32_t max_chunk_vars, Eigen::MatrixXd& geno_mat, Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>& geno_mask);
};

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