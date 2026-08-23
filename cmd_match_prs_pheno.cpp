#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "assoc_utils.h"
#include "qpgen.h"
#include "pheno.h"
#include "qpgen_utils.h"
#include "Eigen/Dense"
#include <cmath>

#include <Eigen/Dense>
#include <cmath>
#include <limits>

struct ShrinkageMinimizer {
    Eigen::MatrixXd MN;          // (M elementwise-product N), p x p
    Eigen::VectorXd d;           // eigenvalues, length p
    Eigen::VectorXd oneMinusD;   // 1 - d_i, length p

    ShrinkageMinimizer(const Eigen::MatrixXd& X,
                       const Eigen::MatrixXd& Y,
                       const Eigen::MatrixXd& U,
                       const Eigen::VectorXd& d_)
        : d(d_),
          oneMinusD(Eigen::VectorXd::Ones(d_.size()) - d_)
    {
        const Eigen::MatrixXd A = X * U;                          // n x p
        const Eigen::MatrixXd B = U.transpose() * Y.transpose();  // p x m
        const Eigen::MatrixXd M = A.transpose() * A;              // p x p (sym)
        const Eigen::MatrixXd N = B * B.transpose();              // p x p (sym)
        MN.noalias() = M.cwiseProduct(N);
    }

    // f(lambda) = ||X U ((1-lambda)D + lambda I)^{-1} U^T Y^T||_F^2
    double squaredNorm(double lambda) const {
        Eigen::VectorXd g = d + lambda * oneMinusD;          // (1-l)d + l
        if ((g.array() <= 0.0).any())
            return std::numeric_limits<double>::infinity();
        Eigen::VectorXd w = g.cwiseInverse();
        double sqnorm = w.dot(MN * w);
        notice("Evaluating squaredNorm at lambda = %g: sqnorm = %g, min(g) = %g, max(g) = %g, min(w) = %g, max(w) = %g", lambda, sqnorm, g.minCoeff(), g.maxCoeff(), w.minCoeff(), w.maxCoeff());
        return sqnorm;
    }

    // Golden-section search on [a, b]. Assumes f is unimodal on the interval
    // (typically true when d_i >= 0 and search range is [0, 1]).
    double findOptimalLambda(double a = 0.0, double b = 1.0,
                             double tol = 1e-9, int maxIter = 200) const {
        const double phi = (std::sqrt(5.0) - 1.0) / 2.0;  // ~0.6180339887
        double x1 = b - phi * (b - a);
        double x2 = a + phi * (b - a);
        double f1 = squaredNorm(x1);
        double f2 = squaredNorm(x2);

        for (int i = 0; i < maxIter && (b - a) > tol; ++i) {
            if (f1 < f2) {
                b  = x2;
                x2 = x1;          f2 = f1;
                x1 = b - phi * (b - a);
                f1 = squaredNorm(x1);
            } else {
                a  = x1;
                x1 = x2;          f1 = f2;
                x2 = a + phi * (b - a);
                f2 = squaredNorm(x2);
            }
        }
        return 0.5 * (a + b);
    }
};

double ledoit_wolf_shrinkage_parameter(const Eigen::MatrixXd& total_cov, const Eigen::MatrixXd& prs_mat, const Eigen::MatrixXd& pheno_mat) {
    int n = prs_mat.rows();
    int p = prs_mat.cols();
    if ( ( p != pheno_mat.cols() ) ) {
        error("PRS and phenotype matrices must have the same dimensions - got %d x %d and %d x %d", (int32_t)prs_mat.rows(), (int32_t)prs_mat.cols(), (int32_t)pheno_mat.rows(), (int32_t)pheno_mat.cols());
    }
    if ( ( p != total_cov.rows() ) || ( p != total_cov.cols() ) ) {
        error("Total covariance matrix must have the same dimensions as the number of traits - got %d x %d", (int32_t)total_cov.rows(), (int32_t)total_cov.cols());
    }

    // perform Eigendecomposition of the total covariance matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver(total_cov);
    if ( eig_solver.info() != Eigen::Success ) {
        error("Failed to perform eigendecomposition of the total covariance matrix");
    }
    Eigen::VectorXd eig_vals = eig_solver.eigenvalues();
    Eigen::MatrixXd eig_vecs = eig_solver.eigenvectors();

    // compute the Ledoit-Wolf shrinkage parameter
    ShrinkageMinimizer minimizer(prs_mat, pheno_mat, eig_vecs, eig_vals);
    double lambda = minimizer.findOptimalLambda(0.0, 1.0);
    double minVal = std::sqrt(minimizer.squaredNorm(lambda));
    notice("Optimal shrinkage parameter (lambda) found: %g with minimum Frobenius norm: %g", lambda, minVal);
    return lambda;
}


int32_t cmd_match_prs_pheno(int32_t argc, char **argv)
{
    std::string prsf;       // PRS files
    std::string phef;
    std::string covf;
    std::string sample_mapf;    // TSV file containing input sample IDs PRS and phenotype files in [PRS_SAMPLE_ID] [PHENO_SAMPLE_ID] format. If they use the same IDs, use only a single column if subsetting samples are needed
    std::string trait_mapf;     // TSV file containing input trait IDs PRS and phenotype files in [PRS_TRAIT_ID] [PHENO_TRAIT_ID] format. If they use the same IDs, use only a single column if subsetting traits are needed
    std::string outf;
    std::string weightf;    // Weights file to use for each phenotype
    std::string prs_format("regenie");
    std::string pheno_format("regenie");
    std::string cov_format("regenie");
    bool rint_after_adj = false;   // Perform rank-based inverse normal transformation after covariate adjustment
    bool use_mahalanobis = false;  // Use Mahalanobis distance for matching
    bool use_ledoit_wolf = false; // Use Ledoit-Wolf shrinkage for covariance estimation when using Mahalanobis distance
    double min_weight = 0.0;  // minimum weight (in r) per trait to set to zero
    double z_lenient_threshold = 1.96;  // Z-score threshold for lenient matching
    double z_diff_threshold = 2.0;      // Z-score difference to declare a clear match
    int32_t n_threads = 1;
    double lambda = 0.0; // regularization parameter for Mahalanobis distance
    double weight_prs_mh = 0.5; // weight for MH distance when combining with weighted correlation

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("prs", &prsf, "Input PRS file")
    LONG_STRING_PARAM("pheno", &phef, "Input phenotype matrix")
    LONG_STRING_PARAM("cov", &covf, "Input covariate matrix (optional)")
    LONG_STRING_PARAM("sample-tsv", &sample_mapf, "TSV file containing input sample IDs PRS and phenotype files in [PRS_SAMPLE_ID] [PHENO_SAMPLE_ID] format. If they use the same IDs, use only a single column if subsetting samples are needed")
    LONG_STRING_PARAM("trait-tsv", &trait_mapf, "TSV file containing input trait IDs PRS and phenotype files in [PRS_TRAIT_ID] [PHENO_TRAIT_ID] format. If they use the same IDs, use only a single column if subsetting traits are needed")
    LONG_STRING_PARAM("weights", &weightf, "Input weights file for each phenotype in [PHENO_ID] [WEIGHT] format")
    LONG_STRING_PARAM("prs-format", &prs_format, "Format of the PRS file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("pheno-format", &pheno_format, "Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("cov-format", &cov_format, "Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")

    LONG_PARAM_GROUP("Analysis options", NULL)
    LONG_PARAM("rint", &rint_after_adj, "Perform rank-based inverse normal transformation after covariate adjustment (default: false)")
    LONG_PARAM("mahalanobis", &use_mahalanobis, "Use Mahalanobis distance for matching (default: false)")
    LONG_PARAM("ledoit-wolf", &use_ledoit_wolf, "Use Ledoit-Wolf shrinkage for covariance estimation when using Mahalanobis distance (default: false)")
    LONG_DOUBLE_PARAM("weight-prs-mh", &weight_prs_mh, "Weight for PRS distance when combining with weighted correlation (default: 0.5)")
    LONG_DOUBLE_PARAM("lambda", &lambda, "Regularization parameter for Mahalanobis distance, between 0 and 1 (default: 0.0)")
    LONG_DOUBLE_PARAM("min-weight", &min_weight, "Minimum weight (in r) per trait to set to zero (default: 0.0)")
    LONG_DOUBLE_PARAM("z-threshold", &z_lenient_threshold, "Z-score threshold for lenient matching (default: 1.96)")
    LONG_DOUBLE_PARAM("z-diff", &z_diff_threshold, "Z-score difference to declare a clear match (default: 2.0)")
    LONG_INT_PARAM("threads", &n_threads, "Number of threads to use (default: 1)")

    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    if ( n_threads > 1 ) {
        Eigen::setNbThreads(n_threads);
        notice("Setting number of threads to %d", n_threads);
    }

    // Load the weight file
    std::map<std::string, double> phe_weights;
    if ( !weightf.empty() ) {
        notice("Loading phenotype weights from %s", weightf.c_str());
        tsv_reader tr_weight(weightf.c_str());
        while( tr_weight.read_line() ) {
            if ( tr_weight.nfields < 2 ) {
                error("Invalid format weights file %s in line %zu starting with %s. Must containing at least 2 fields", weightf.c_str(), (int32_t)weightf.size() + 1, tr_weight.str_field_at(0) );
            }
            std::string phe_id = tr_weight.str_field_at(0);
            double weight = tr_weight.double_field_at(1);
            phe_weights[phe_id] = weight;
        }
        if ( phe_weights.size() == 0 ) {
            error("No valid phenotype weights found in file %s", weightf.c_str());
        }
    }

    notice("Loading phenotype matrix from %s", phef.c_str());

    // load the phenotype matrix
    PhenoMatrix pheno_matrix;
    if ( !pheno_matrix.load_pheno_matrix(phef.c_str(), pheno_format.c_str()) ) {
        error("Failed to load the phenotype matrix from file %s", phef.c_str());
    }

    notice("Loaded phenotype matrix with %d samples and %d phenotypes from %s", (int32_t)pheno_matrix.samp_ids.size(), (int32_t)pheno_matrix.pheno_ids.size(), phef.c_str());

    notice("Loading PRS matrix from %s", prsf.c_str());
    // load the PRS matrix
    PhenoMatrix prs_matrix;
    if ( !prs_matrix.load_pheno_matrix(prsf.c_str(), prs_format.c_str()) ) {
        error("Failed to load the PRS matrix from file %s", prsf.c_str());
    }

    notice("Loaded PRS matrix with %d samples and %d traits from %s", (int32_t)prs_matrix.samp_ids.size(), (int32_t)prs_matrix.pheno_ids.size(), prsf.c_str());

    // load the covariance matrix
    PhenoMatrix cov_matrix;
    if ( !covf.empty() ) {
        notice("Loading covariate matrix from %s", covf.c_str());
        if ( !cov_matrix.load_pheno_matrix(covf.c_str(), cov_format.c_str()) ) {
            error("Failed to load the covariate matrix from file %s", covf.c_str());
        }
        notice("Loaded covariate matrix with %d samples and %d covariates from %s", (int32_t)cov_matrix.samp_ids.size(), (int32_t)cov_matrix.pheno_ids.size(), covf.c_str());
    }

    // subset to overlapping phenotypes
    std::vector<int32_t> phe_pheno_indices;
    std::vector<int32_t> prs_pheno_indices;
    if ( trait_mapf.empty() ) { // identify overlapping phenotypes based on names
        std::vector<std::string> overlapping_phenos;
        identify_overlapping_ids(prs_matrix.pheno_ids, pheno_matrix.pheno_ids, overlapping_phenos);
        notice("Found %d overlapping phenotypes between PRS and phenotype matrices based on phenotype IDs", (int32_t)overlapping_phenos.size());
        if ( overlapping_phenos.size() == 0 ) {
            error("No overlapping phenotypes found between PRS and phenotype matrices. Please check if the phenotype IDs are consistent across the files or provide a trait ID mapping file using --trait-tsv option.");
        }
        pheno_matrix.subset_pheno_ids(overlapping_phenos);
        prs_matrix.subset_pheno_ids(overlapping_phenos);
    }
    else {
        notice("Loading trait ID mapping between PRS and phenotype files from %s", trait_mapf.c_str());
        tsv_reader tr_trait_map(trait_mapf.c_str());
        while( tr_trait_map.read_line() ) {
            if ( tr_trait_map.nfields > 2 ) {
                error("Invalid format trait ID mapping file %s in line %zu starting with %s. Must containing at least 2 fields", trait_mapf.c_str(), (int32_t)trait_mapf.size() + 1, tr_trait_map.str_field_at(0) );
            }
            std::string prs_id = tr_trait_map.str_field_at(0);
            std::string phe_id = tr_trait_map.str_field_at(tr_trait_map.nfields == 1 ? 0 : 1);
            std::map<std::string, int32_t>::const_iterator it_prs = prs_matrix.pheno_id2idx.find(prs_id);
            std::map<std::string, int32_t>::const_iterator it_phe = pheno_matrix.pheno_id2idx.find(phe_id);
            if ( it_prs == prs_matrix.pheno_id2idx.end() ) {
                error("PRS trait ID %s in mapping file %s not found in PRS matrix", prs_id.c_str(), trait_mapf.c_str());
            }
            if ( it_phe == pheno_matrix.pheno_id2idx.end() ) {
                error("Phenotype trait ID %s in mapping file %s not found in phenotype matrix", phe_id.c_str(), trait_mapf.c_str());
            }
            prs_pheno_indices.push_back( it_prs->second );
            phe_pheno_indices.push_back( it_phe->second );
        }
        pheno_matrix.subset_pheno_indices( phe_pheno_indices );
        prs_matrix.subset_pheno_indices( prs_pheno_indices );
    }

    // after the subsetting, the number of phenotypes in both matrices should be the same
    if ( pheno_matrix.pheno_ids.size() != prs_matrix.pheno_ids.size() ) {
        error("Number of phenotypes after subsetting do not match between phenotype and PRS matrices (%d vs %d)", (int32_t)pheno_matrix.pheno_ids.size(), (int32_t)prs_matrix.pheno_ids.size());
    }
    notice("Subsetted to %d overlapping phenotypes between PRS and phenotype matrices", (int32_t)pheno_matrix.pheno_ids.size());

    // if covariates are provided, subset to overlapping samples
    if ( !covf.empty() ) {
        std::vector<std::string> overlapping_samples;
        identify_overlapping_ids(pheno_matrix.samp_ids, cov_matrix.samp_ids, overlapping_samples);
        notice("Found %d overlapping samples between phenotype and covariate matrices", (int32_t)overlapping_samples.size());
        if ( overlapping_samples.size() == 0 ) {
            error("No overlapping samples found between phenotype and covariate matrices. Please check if the sample IDs are consistent across the files.");
        }
        pheno_matrix.subset_sample_ids(overlapping_samples);
        cov_matrix.subset_sample_ids(overlapping_samples);
        notice("Subsetted phenotype and covariate matrices to %d overlapping samples", (int32_t)overlapping_samples.size());

        notice("Adjusting phenotypes by covariates using linear regression");
        if ( pheno_matrix.has_missing || cov_matrix.has_missing ) {
            error("Covariate adjustment is currently only supported for phenotype and covariate matrices without missing values");
        }
        else {
            pheno_matrix.pheno_mat = pheno_adj_cov_nxt_without_missing(pheno_matrix.pheno_mat, cov_matrix.pheno_mat);
        }
    }

    if ( rint_after_adj ) {
        notice("Performing rank-based inverse normal transformation for all phenotypes after covariate adjustment");
        if ( !pheno_matrix.has_missing ) {
            pheno_matrix.pheno_mat = rint_matrix_without_missing(pheno_matrix.pheno_mat);
        }
        else {
            error("Rank-based inverse normal transformation adjustment is currently only supported for phenotype matrices with missing values");
        }
    }

    // construct matching ID maps
    std::vector<int32_t> matching_prs_samp_indices;
    std::vector<int32_t> matching_pheno_samp_indices;
    std::map<int32_t, int32_t> samp_idx_pheno2prs;
    bool has_sample_map = !sample_mapf.empty();
    bool has_trait_map = !trait_mapf.empty();

    // identify PRS samples to use
    if ( has_sample_map ) {
        notice("Loading sample ID mapping between PRS and phenotype files from %s", sample_mapf.c_str());
        tsv_reader tr_sample_map(sample_mapf.c_str());
        while( tr_sample_map.read_line() ) {
            if ( tr_sample_map.nfields > 2 ) {
                error("Invalid format sample ID mapping file %s in line %zu starting with %s. Must containing at least 2 fields", sample_mapf.c_str(), (int32_t)sample_mapf.size() + 1, tr_sample_map.str_field_at(0) );
            }
            std::string prs_id = tr_sample_map.str_field_at(0);
            std::string phe_id = tr_sample_map.str_field_at(tr_sample_map.nfields == 1 ? 0 : 1);
            std::map<std::string, int32_t>::const_iterator it_prs = prs_matrix.samp_id2idx.find(prs_id);
            std::map<std::string, int32_t>::const_iterator it_phe = pheno_matrix.samp_id2idx.find(phe_id);
            if ( it_prs == prs_matrix.samp_id2idx.end() ) {
                error("PRS trait ID %s in mapping file %s not found in PRS matrix", prs_id.c_str(), sample_mapf.c_str());
            }
            if ( it_phe == pheno_matrix.samp_id2idx.end() ) {
                error("Phenotype trait ID %s in mapping file %s not found in phenotype matrix", phe_id.c_str(), sample_mapf.c_str());
            }
            matching_prs_samp_indices.push_back( it_prs->second );
            matching_pheno_samp_indices.push_back( it_phe->second );
            samp_idx_pheno2prs[it_phe->second] = it_prs->second;
        }
        notice("Found %d overlapping samples between PRS and phenotype matrices based on sample ID mapping file", (int32_t)matching_pheno_samp_indices.size());
    }
    else { // identify overlapping samples based on names
        for(int32_t i=0; i < pheno_matrix.samp_ids.size(); ++i ) {
            const std::string& samp_id = pheno_matrix.samp_ids[i];
            std::map<std::string, int32_t>::const_iterator it = prs_matrix.samp_id2idx.find( samp_id );
            if ( it != prs_matrix.samp_id2idx.end() ) {
                matching_pheno_samp_indices.push_back( i );
                matching_prs_samp_indices.push_back( it->second );
                samp_idx_pheno2prs[i] = it->second;
            }
        }
        notice("Found %d overlapping samples between PRS and phenotype matrices based on sample IDs", (int32_t)matching_pheno_samp_indices.size());
    }

    // standardize each trait
    notice("Standardizing PRS and phenotype matrices");
    standardize_matrix_columns_inplace(prs_matrix.pheno_mat);
    standardize_matrix_columns_inplace(pheno_matrix.pheno_mat);

    // for PRS X, and phenotype Y, compute column-wise correlation between X and Y
    Eigen::VectorXd weights( pheno_matrix.pheno_ids.size() );
    if ( weightf.empty() ) {
        notice("Computing weights for each phenotype");
        //Eigen::VectorXd weights = columnwise_dot(prs_matrix.pheno_mat, pheno_matrix.pheno_mat) / (double)n_overlapping_samples;
        for(int32_t i=0; i < pheno_matrix.pheno_ids.size(); ++i) {
            double r = 0;
            for(int32_t j=0; j < matching_pheno_samp_indices.size(); ++j) {
                r += prs_matrix.pheno_mat( matching_prs_samp_indices[j], i ) * pheno_matrix.pheno_mat( matching_pheno_samp_indices[j], i );
            }
            weights[i] = r / (double)matching_pheno_samp_indices.size();
        }

        htsFile* wf = hts_open((outf + ".weights.tsv.gz").c_str(), "wz");
        if ( wf == NULL ) {
            error("Cannot open output file %s.weights.tsv for writing", outf.c_str());
        }
        hprintf(wf, "Trait\tWeight\n");
        for ( int32_t i = 0; i < weights.size(); ++i) {
            hprintf(wf, "%s\t%.6f\n", pheno_matrix.pheno_ids[i].c_str(), weights[i]);
        }
        hts_close(wf);
    }
    else {
        notice("Using provided weights for each phenotype from %s", weightf.c_str());
        for(int32_t i=0; i < pheno_matrix.pheno_ids.size(); ++i) {
            std::map<std::string, double>::const_iterator it = phe_weights.find( pheno_matrix.pheno_ids[i] );
            if ( it == phe_weights.end() ) {
                notice("No weight found for phenotype %s in weights file %s, setting weight to 0", pheno_matrix.pheno_ids[i].c_str(), weightf.c_str());
                weights[i] = 0.0;
            }
            else {
                weights[i] = it->second;
            }
        }
    }

    // set weights with abs value < min_weight to zero
    int32_t n_pass_weights = 0;
    for ( int32_t i = 0; i < weights.size(); ++i ) {
        if ( weights[i] < min_weight ) {
            weights[i] = 0.0;
        }
        else {
            n_pass_weights++;
        }
    }
    notice("%d / %zu phenotypes passed the minimum weight threshold of %.4f", n_pass_weights, weights.size(), min_weight);

    // calculate all pair weighted correlations [n_prs x n_pheno] matrix
    Eigen::MatrixXd all_pair_wcor;
    if ( use_mahalanobis ) {
        notice("Using Mahalanobis distance for matching to account for correlation between traits (may take much longer time)");
        notice("Computing PRS covariance matrix");
        Eigen::MatrixXd prs_cov = prs_matrix.pheno_mat.transpose() * prs_matrix.pheno_mat / (double)prs_matrix.pheno_mat.rows();
        notice("Computing phenotype covariance matrix");
        Eigen::MatrixXd pheno_cov = pheno_matrix.pheno_mat.transpose() * pheno_matrix.pheno_mat / (double)pheno_matrix.pheno_mat.rows();
        notice("Comptuing total covariance matrix");
        Eigen::MatrixXd total_cov = weight_prs_mh * prs_cov + (1.0 - weight_prs_mh) * pheno_cov + 1e-8 * Eigen::MatrixXd::Identity( prs_cov.rows(), prs_cov.cols() );
        if ( lambda < 0.0 || lambda > 1.0 ) {
            error("Invalid value for lambda: %.4f. Must be between 0 and 1", lambda);
        }
        if ( use_ledoit_wolf ) {
            if ( lambda > 0 ) {
                error("Cannot use Ledoit-Wolf shrinkage when lambda is greater than 0. Please set lambda to 0 to use Ledoit-Wolf shrinkage or set lambda to a value between 0 and 1 to use regularization without Ledoit-Wolf shrinkage");
            }
            else {
                lambda = ledoit_wolf_shrinkage_parameter(total_cov, prs_matrix.pheno_mat, pheno_matrix.pheno_mat);
                notice("Using Ledoit-Wolf shrinkage with lambda = %.4f", lambda);
            }
        }
        if ( lambda > 0 ) {
            notice("Adding regularization parameter %.4f to total covariance matrix", lambda);
            total_cov = (1-lambda) * total_cov + lambda * Eigen::MatrixXd::Identity( prs_cov.rows(), prs_cov.cols() );
        }
        notice("Inverting total covariance matrix");
        Eigen::MatrixXd total_cov_inv = total_cov.inverse();

        notice("Constructing weighted Mahalanobis matrix to multiply");
        Eigen::MatrixXd W_sqrt = weights.cwiseSqrt().asDiagonal();
        Eigen::MatrixXd M = W_sqrt * total_cov_inv * W_sqrt;

        notice("Computing numerator for all pair weighted correlations between PRS traits and phenotypes");
        Eigen::MatrixXd numerator = prs_matrix.pheno_mat * M * pheno_matrix.pheno_mat.transpose();

        notice("Computing norms for PRS and phenotype matrices");
        Eigen::VectorXd prs_norms = (prs_matrix.pheno_mat * M * prs_matrix.pheno_mat.transpose()).diagonal().cwiseSqrt();
        Eigen::VectorXd pheno_norms = (pheno_matrix.pheno_mat * M * pheno_matrix.pheno_mat.transpose()).diagonal().cwiseSqrt();

        // 4. Compute Correlation Matrix
        // Divide numerator(i, j) by (prs_norm(i) * pheno_norm(j))
        // Using broadcasting or a loop (Eigen broadcasting can be tricky, loop is safe)
        all_pair_wcor.resize(numerator.rows(), numerator.cols());
        notice("Computing all pair weighted correlations between PRS traits and phenotypes");
        for (int i = 0; i < numerator.rows(); ++i) {
            for (int j = 0; j < numerator.cols(); ++j) {
                double denom = prs_norms(i) * pheno_norms(j);
                all_pair_wcor(i, j) = numerator(i, j) / (denom + 1e-100);
            }
        }
        //notice("Computing all pair Mahalanobis distances between PRS traits and phenotypes");
        //all_pair_wcor = ( prs_matrix.pheno_mat * ( ( W_sqrt * total_cov_inv * W_sqrt ) * pheno_matrix.pheno_mat.transpose() ) ) / ( weights.array().abs().sum() + 1e-100 );

        // notice("Computing difference between PRS and phenotype matrices");
        // Eigen::MatrixXd diff_mat = prs_matrix.pheno_mat - pheno_matrix.pheno_mat;

        // notice("Computing all pair Mahalanobis distances between PRS traits and phenotypes");
        // all_pair_wcor =  -1 * ( diff_mat * ( ( total_cov_inv * weights.asDiagonal() ) * diff_mat.transpose() ) );
    }
    else {
        notice("Computing all pair weighted correlations between PRS traits and phenotypes assuming independence");
        all_pair_wcor = prs_matrix.pheno_mat * ( weights.asDiagonal() * pheno_matrix.pheno_mat.transpose() ) / ( weights.array().abs().sum() + 1e-100 );
    }

    notice("Standardizing all pair weighted correlation matrix");
    // copy the weighted correlation matrix
    Eigen::MatrixXd all_pair_z = all_pair_wcor;
    // Apply column-wise standardization
    standardize_matrix_columns_inplace(all_pair_z);

    // open the output file gz or plain based on the extension
    htsFile* wf1 = hts_open((outf + ".match.assigned.tsv.gz").c_str(), "wz");
    if ( wf1 == NULL ) {
        error("Cannot open output file %s.match.assigned.tsv.gz for writing", outf.c_str());
    }
    htsFile* wf2 = hts_open((outf + ".match.all.tsv.gz").c_str(), "wz");
    if ( wf2 == NULL ) {
        error("Cannot open output file %s.match.all.tsv.gz for writing", outf.c_str());
    }
 
    // write the header line
    hprintf(wf1, "ID.Pheno\tMatchStatus\tID.self\tZ.self\tCOR.self\tRank.self\n");
    hprintf(wf2, "ID.Pheno\tMatchStatus\tID.self\tZ.self\tCOR.self\tRank.self\tID.1st\tZ.1st\tCOR.1st\tID.2nd\tZ.2nd\tCOR.2nd\tID.3rd\tZ.3rd\tCOR.3rd\tID.4th\tZ.4th\tCOR.4th\tID.5th\tZ.5th\tCOR.5th\n");
    for ( int32_t i = 0; i < pheno_matrix.samp_ids.size(); ++i ) {
        // find the best and second best matches
        //notice("foo %d", i); 
        int32_t self_idx = -1;
        if ( samp_idx_pheno2prs.find(i) != samp_idx_pheno2prs.end() ) {
            self_idx = samp_idx_pheno2prs[i];
        }
        int32_t self_rank = 1;
        double z_self = self_idx >=0 ? all_pair_z(self_idx, i) : -9999.0;
        double cor_self = self_idx >=0 ? all_pair_wcor(self_idx, i) : -9999.0;
        double idx_top[5] = { -1, -1, -1, -1, -1 };
        double z_top[5] = { -9999.0, -9999.0, -9999.0, -9999.0, -9999.0 };
        double cor_top[5] = { -9999.0, -9999.0, -9999.0, -9999.0, -9999.0 };
        for( int32_t j = 0; j < prs_matrix.samp_ids.size(); ++j ) {
            double z = all_pair_z(j, i);
            double cor = all_pair_wcor(j, i);
            if ( self_idx >= 0 && j != self_idx && z > z_self ) {
                self_rank++;
            }
            // check if this is in the top 5
            for ( int32_t k = 0; k < 5; ++k ) {
                if ( z > z_top[k] ) {
                    // shift down
                    for ( int32_t l = 4; l > k; --l ) {
                        z_top[l] = z_top[l-1];
                        cor_top[l] = cor_top[l-1];
                        idx_top[l] = idx_top[l-1];
                    }
                    z_top[k] = z;
                    cor_top[k] = cor;
                    idx_top[k] = j;
                    break;
                }
            }
        }
        if ( self_idx >= 0 ) {
            const char* match_status = ( self_rank == 1 ? "BEST_MATCH" : ( z_self < z_lenient_threshold ? "NO_MATCH" : "LENIENT_MATCH" ) );
            hprintf(wf1, "%s\t%s\t%s\t%.6f\t%.6f\t%d\n",
                pheno_matrix.samp_ids[i].c_str(),
                match_status,
                prs_matrix.samp_ids[self_idx].c_str(),
                z_self, cor_self,
                self_rank);
        }
        const char* match_status = "UNCLEAR";
        if ( self_rank == 1 ) {
            match_status = "SELF_BEST";
        }
        else if ( z_top[0] > z_top[1] + z_diff_threshold ) {
            match_status = "SINGLE_NEW_BEST";
        }
        else {
            if ( z_self > z_lenient_threshold ) {
                match_status = "SELF_LENIENT";
            }
            for(int32_t k=2; k < 5; ++k ) {
                if ( z_top[k-1] > z_top[k] + z_diff_threshold ) {
                    match_status = "MULTI_NEW_BEST";
                    break;
                }
            }
        }
        if ( self_idx >= 0 ) {
            hprintf(wf2, "%s\t%s\t%s\t%.6f\t%.6f\t%d",
                pheno_matrix.samp_ids[i].c_str(),
                match_status,
                prs_matrix.samp_ids[self_idx].c_str(),
                z_self, cor_self,
                self_rank);
        }
        else {
            hprintf(wf2, "%s\t%s\tNA\tNA\tNA\tNA",
                pheno_matrix.samp_ids[i].c_str(),
                match_status);
        }
        for ( int32_t k = 0; k < 5; ++k ) {
            if ( idx_top[k] >=0 ) {
                hprintf(wf2, "\t%s\t%.6f\t%.6f",
                    prs_matrix.samp_ids[ (int32_t)idx_top[k] ].c_str(),
                    z_top[k],
                    cor_top[k]);
            }
            else {
                hprintf(wf2, "\tNA\tNA\tNA");
            }
        }
        hprintf(wf2, "\n");
    }
    hts_close(wf1); // close the output file
    hts_close(wf2); // close the output file

 
    notice("Analysis finished");
    return 0;
}