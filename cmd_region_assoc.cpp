#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "assoc_utils.h"
#include "qpgen_utils.h"
#include "qpgen.h"
#include "pheno.h"
#include "Eigen/Dense"
#include <cmath>

int32_t cmd_region_assoc(int32_t argc, char **argv)
{
    std::string pgenf;      // PLINK2 genotype file
    std::string psamf;      // PLINK2 sample file
    std::string pivarf;     // PLINK2 indexed pvar file, compatible with PgenIdxReader
    std::string pgenlistf;  // CHROM BEG END PGEN PSAM PIVAR containing the list of region-specific BED files
    std::string phef;
    std::string covf;
    std::string outf;
    std::string samplef;

    // information about the trait and variants to be tested
    std::string traits;   // assume that a single trait is being tested
    std::string traitf;
    std::string region;
    double min_af = 0.0;
    double max_af = 1.0;
    double min_ac = 0.0;
    double max_ac = 1e9;

    int32_t jump_thres_bp = 1000000; 
    int32_t max_chunk_vars = 1000; // Maximum number of variants to store at once in memory
    int32_t offset_pheno = 1; 
    int32_t offset_cov = 1;
    int32_t icol_pivar_idx = 9;
    int32_t icol_pheno_id = 1;
    int32_t icol_cov_id = 1;
    std::string pheno_format("regenie");
    std::string cov_format("regenie");
    std::string colname_pheno_sample("pheno"); // column name for the sample IDs in the phenotype file
    std::string colname_geno_sample("geno"); // column name for the sample IDs in the genotype file 
    bool rint_before_adj = false; // Perform rank-based inverse normal transformation before covariate adjustment
    bool rint_after_adj = false; // Perform rank-based inverse normal transformation after covariate adjustment

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("pgen", &pgenf, "Input PLINK2 genotype file")
    LONG_STRING_PARAM("psam", &psamf, "Input PLINK2 sample file")
    LONG_STRING_PARAM("pivar", &pivarf, "Input PLINK2 index pvar file (bgzipped and tabix)")
    LONG_STRING_PARAM("pgen-list", &pgenlistf, "Input file containing CHROM BEG END PGEN PSAM PIVAR")
    LONG_STRING_PARAM("pheno", &phef, "Input phenotype matrix")
    LONG_STRING_PARAM("sample", &samplef, "Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files")
    LONG_STRING_PARAM("cov", &covf, "Input covariate matrix (optional)")
    LONG_STRING_PARAM("pheno-format", &pheno_format, "Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("cov-format", &cov_format, "Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("traits", &traits, "Trait IDs (comma-separated) to be tested (required)")
    LONG_STRING_PARAM("traitf", &traitf, "Input file containing trait IDs to be tested (one per line)")
    LONG_STRING_PARAM("region", &region, "Region string in the format CHROM:BEG-END (required)")
    LONG_DOUBLE_PARAM("min-af", &min_af, "Minimum allele frequency for variants to be tested (default: 0.0)")
    LONG_DOUBLE_PARAM("max-af", &max_af, "Maximum allele frequency for variants to be tested (default: 1.0)")
    LONG_DOUBLE_PARAM("min-ac", &min_ac, "Minimum allele count for variants to be tested (default: 0.0)")
    LONG_DOUBLE_PARAM("max-ac", &max_ac, "Maximum allele count for variants to be tested (default: 1e9)")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")

    LONG_PARAM_GROUP("Auxiliary options", NULL)
    LONG_INT_PARAM("jump-thres-bp", &jump_thres_bp, "Jump threshold in base pairs for the variant index (default: 1000000)")
    LONG_INT_PARAM("max-chunk-vars", &max_chunk_vars, "Maximum number of variants to store at once in memory (default: 1000)")
    LONG_INT_PARAM("offset-pheno", &offset_pheno, "The number of index columns (i.e. not containing values) in the phenotype files (default: 1)")
    LONG_INT_PARAM("offset-cov", &offset_cov, "The number of index columns (i.e. not containing values) in the covariate columns (default: 1)")
    LONG_INT_PARAM("icol-pivar-idx", &icol_pivar_idx, "1-based column index for the variant ID in the pvar file (default: 9)")
    LONG_INT_PARAM("icol-pheno-id", &icol_pheno_id, "1-based column index for the phenotype ID in the phenotype file (default: 1)")
    LONG_INT_PARAM("icol-cov-id", &icol_cov_id, "1-based column index for the covariate ID in the phenotype file (default: 1)")
    LONG_STRING_PARAM("colname-pheno-sample", &colname_pheno_sample, "When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')")
    LONG_STRING_PARAM("colname-geno-sample", &colname_geno_sample, "When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'geno')")
    LONG_PARAM("rint-before-adj", &rint_before_adj, "Perform rank-based inverse normal transformation before covariate adjustment (default: false)")
    LONG_PARAM("rint-after-adj", &rint_after_adj, "Perform rank-based inverse normal transformation after covariate adjustment (default: false)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // check required arguments
    if ( phef.empty() ) {
        error("Phenotype file (--pheno) is required");
    }
    if ( region.empty() ) {
        error("Region (--region) is required");
    }
    
    // set up the association input object and configure the parameters
    ind_assoc_input input;
    input.set_jump_thres_bp(jump_thres_bp);
    input.set_icol_pivar_idx(icol_pivar_idx - 1); // convert to 0-based index
    input.set_rint_before_adj(rint_before_adj);
    input.set_rint_after_adj(rint_after_adj);
    input.set_minmax_af(min_af, max_af);
    input.set_minmax_ac(min_ac, max_ac);
    if ( !samplef.empty() ) {
        input.set_subset_sample_file(samplef.c_str());
    }
    if ( !traitf.empty() ) {
        input.set_subset_pheno_file(traitf.c_str(), icol_pheno_id);
        if ( !traits.empty() ) {
            error("Cannot provide both --traits and --traitf");
        }
    }
    else if ( !traits.empty() ) {
        // split the trait string by comma and set the subset phenotype IDs
        std::vector<std::string> subset_pheno_ids;
        split(subset_pheno_ids, ",", traits);
        input.set_subset_pheno_ids(subset_pheno_ids);
    }
    else {
        error("Either --traits or --traitf must be provided");
    }

    notice("Loading genotype data from pgen files with chromosome %s, position %d to %d, and maximum chunk size of %d variants", region.c_str(), 0, 0, max_chunk_vars);
    if ( !pgenlistf.empty() ) { // list is provided
        input.process_pgenlist(pgenlistf.c_str(), phef.c_str(), pheno_format.c_str(), covf.c_str(), cov_format.c_str());
    }
    else if ( !pgenf.empty() && !pivarf.empty() && !psamf.empty() ) { // single pgen file is provided
        input.process_single_pgen(pgenf.c_str(), pivarf.c_str(), psamf.c_str(), phef.c_str(), pheno_format.c_str(), covf.c_str(), cov_format.c_str());
    }
    else {
        error("Either --pgen-list or --pgen, --pivar, and --psam must be provided");
    }    

    cbe_t region_cbe(region.c_str());
    // read chunk of genotype based on the region
    Eigen::MatrixXd geno_mat;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> geno_mask;
    notice("icol_pivar_idx: %d", input.mpr.get_icol_pivar_idx());
    bool first = true;
    bool any_loaded = false;
    int64_t total_vars = 0;
    // load_genotype_chunk returns true (and fills geno_mat) for every chunk with >=1 variant,
    // and false only once the region is fully streamed, so each loaded chunk is processed exactly once.
    while ( input.load_genotype_chunk(first ? region_cbe.chrom.c_str() : NULL, region_cbe.beg1, region_cbe.end0, max_chunk_vars, geno_mat, geno_mask) ) {
        first = false;
        any_loaded = true;
        total_vars += geno_mat.cols();

        // perform rectangular association analysis
        // std::vector<std::vector<slr_sumstat_t> > rect_results;
        // notice("Performing rectangular association analysis for %d phenotypes and %d variants after skipping %d variants", (int32_t)pheno_matrix.pheno_ids.size(), new_chunk_size, n_skipped);
        // if ( !simple_rect_regression_without_missing(
        //         input.pheno_matrix.pheno_mat,
        //         geno_mat,
        //         rect_results) ) {
        //     error("Failed to perform rectangular association analysis for chunk %d", i / max_chunk_vars + 1);
        // }
        // // write the results
        // for(int32_t j=0; j < new_chunk_size; ++j) {
        //     cpra_t cpra(chunk_cpras[j].c_str());
        //     const var_cnt_t& vcnt = chunk_var_cnts[j];
        //     double a1freq = (double)(vcnt.ac) / (double)(vcnt.an);
        //     hprintf(wf, "%s\t%d\t%s\t%s\t%s\t%.6g\t%d\t%d\t%d\t%d\tLinear",
        //         cpra.chrom.c_str(),
        //         cpra.pos,
        //         chunk_cpras[j].c_str(),
        //         cpra.ref.c_str(),
        //         cpra.alts.c_str(),
        //         a1freq,
        //         vcnt.an,
        //         vcnt.gcs[0],
        //         vcnt.gcs[1],
        //         vcnt.gcs[2]);
        //     const std::vector<slr_sumstat_t>& var_results = rect_results[j];
        //     for(int32_t k=0; k < (int32_t)pheno_matrix.pheno_ids.size(); ++k) {
        //         const slr_sumstat_t& ss = var_results[k];
        //         hprintf(wf, "\t%.6g\t%.6g\t%.6g\t%.6g",
        //             ss.beta,
        //             ss.se,
        //             ss.tstat,
        //             ss.log10p);
        //     }
        //     hprintf(wf, "\n");
        // }

        notice("Performing association mapping for %d variants in the region %s:%d-%d", (int32_t)geno_mat.cols(), region_cbe.chrom.c_str(), region_cbe.beg1, region_cbe.end0);
    }
    if ( !any_loaded ) {
        notice("No variants found in the specified region %s:%d-%d", region_cbe.chrom.c_str(), region_cbe.beg1, region_cbe.end0);
    }
    notice("Analysis finished. Total variants processed: %lld", (long long)total_vars);
    return 0;
}