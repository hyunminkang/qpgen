#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "assoc_utils.h"
#include "susie_utils.h"
#include "qpgen_utils.h"
#include "qpgen.h"
#include "pheno.h"
#include "Eigen/Dense"
#include <cmath>
#include <string>
#include <algorithm>

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
    //int32_t max_chunk_vars = 1000; // Maximum number of variants to store at once in memory
    int32_t max_allowed_vars = 1000000; // Maximum number of variants to store at once in memory
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

    // SuSiE fine-mapping options
    bool run_susie = false;
    int32_t susie_L = 10;
    int32_t susie_max_iter = 100;
    double susie_coverage = 0.95;
    double susie_min_abs_corr = 0.5;
    double susie_tol = 1e-3;
    bool susie_no_standardize = false;
    bool output_lbf = false; // also write per-variant log Bayes factors
    std::string unmappable_effects = "none"; // "none" (standard SuSiE) or "inf" (SuSiE-inf)
    std::string ash_fix_pi_str;              // reduction-test: comma-separated pi values
    std::string ash_fix_sa2_str;             // reduction-test: comma-separated sa2 values

    // suffix for the output files
    std::string assoc_suffix = ".assoc.tsv.gz";
    std::string susie_cs_suffix = ".susie.cs.tsv.gz";
    std::string susie_lbf_suffix = ".susie.lbf.tsv.gz";

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
    LONG_STRING_PARAM("assoc-suffix", &assoc_suffix, "Suffix for the association output file (default: '.assoc.tsv.gz')")
    LONG_STRING_PARAM("susie-cs-suffix", &susie_cs_suffix, "Suffix for the SuSiE credible set output file (default: '.susie.cs.tsv.gz')")
    LONG_STRING_PARAM("susie-lbf-suffix", &susie_lbf_suffix, "Suffix for the SuSiE log Bayes factor output file (default: '.susie.lbf.tsv.gz')")

    LONG_PARAM_GROUP("Auxiliary options", NULL)
    LONG_INT_PARAM("jump-thres-bp", &jump_thres_bp, "Jump threshold in base pairs for the variant index (default: 1000000)")
    LONG_INT_PARAM("max-allowed-vars", &max_allowed_vars, "Maximum number of allowed variants to store at once in memory (default: 1000000)")
    LONG_INT_PARAM("offset-pheno", &offset_pheno, "The number of index columns (i.e. not containing values) in the phenotype files (default: 1)")
    LONG_INT_PARAM("offset-cov", &offset_cov, "The number of index columns (i.e. not containing values) in the covariate columns (default: 1)")
    LONG_INT_PARAM("icol-pivar-idx", &icol_pivar_idx, "1-based column index for the variant ID in the pvar file (default: 9)")
    LONG_INT_PARAM("icol-pheno-id", &icol_pheno_id, "1-based column index for the phenotype ID in the phenotype file (default: 1)")
    LONG_INT_PARAM("icol-cov-id", &icol_cov_id, "1-based column index for the covariate ID in the phenotype file (default: 1)")
    LONG_STRING_PARAM("colname-pheno-sample", &colname_pheno_sample, "When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'pheno')")
    LONG_STRING_PARAM("colname-geno-sample", &colname_geno_sample, "When --sample is provided, the column name for the sample IDs in the phenotype file (default: 'geno')")
    LONG_PARAM("rint-before-adj", &rint_before_adj, "Perform rank-based inverse normal transformation before covariate adjustment (default: false)")
    LONG_PARAM("rint-after-adj", &rint_after_adj, "Perform rank-based inverse normal transformation after covariate adjustment (default: false)")

    LONG_PARAM_GROUP("SuSiE fine-mapping options", NULL)
    LONG_PARAM("susie", &run_susie, "Run SuSiE fine-mapping for each tested phenotype in the region")
    LONG_INT_PARAM("susie-L", &susie_L, "Maximum number of causal single effects (default: 10)")
    LONG_INT_PARAM("susie-max-iter", &susie_max_iter, "Maximum number of IBSS iterations (default: 100)")
    LONG_DOUBLE_PARAM("susie-coverage", &susie_coverage, "Target coverage of credible sets (default: 0.95)")
    LONG_DOUBLE_PARAM("susie-min-abs-corr", &susie_min_abs_corr, "Minimum absolute correlation (purity) required to report a credible set (default: 0.5)")
    LONG_DOUBLE_PARAM("susie-tol", &susie_tol, "Convergence tolerance for the SuSiE objective (default: 1e-3)")
    LONG_PARAM("susie-no-standardize", &susie_no_standardize, "Do not standardize genotype columns to unit variance before SuSiE")
    LONG_STRING_PARAM("unmappable-effects", &unmappable_effects, "Unmappable-effects model for SuSiE: 'none' (standard), 'inf' (SuSiE-inf, adds an infinitesimal effect), or 'ash' (SuSiE-ash, scale-mixture prior). Matches run_susie_v1.r --method (default: none)")
    LONG_STRING_PARAM("ash-fix-pi", &ash_fix_pi_str, "Reduction test: comma-separated pi vector to hold ash mixture weights fixed (skips EM). Length K. Used with --unmappable-effects ash to prove ash reduces to none (pi=1,0,...,0) or inf (pi=0,...,0,1).")
    LONG_STRING_PARAM("ash-fix-sa2", &ash_fix_sa2_str, "Reduction test: comma-separated sa2 grid to hold ash prior-variance grid fixed. Length must match --ash-fix-pi.")
    LONG_PARAM("output-lbf", &output_lbf, "Also write per-variant log Bayes factors (one column per single effect) when running SuSiE")
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

    notice("Loading genotype data from pgen files with chromosome %s, position %d to %d, and maximum chunk size of %d variants", region.c_str(), 0, 0, max_allowed_vars);
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
    // bool first = true;
    // bool any_loaded = false;
    // int64_t total_vars = 0;

    // read all genotypes as one single chunk
    int32_t n_loaded_vars = input.load_genotype_chunk(region_cbe.chrom.c_str(), region_cbe.beg1, region_cbe.end0, max_allowed_vars);
    if ( n_loaded_vars >= max_allowed_vars ) {
        error("Maximum number of allowed variants reached. Please increase the --max-chunk-vars parameter to a larger value (current value: %d)");
    }

    // perform rectangular association analysis
    std::vector<std::vector<slr_sumstat_t> > rect_results;
    notice("Performing rectangular association analysis...");
    if ( !simple_rect_regression_without_missing(
            input.pheno_matrix.pheno_mat,
            input.geno_chunk.geno_mat,
            rect_results) ) {
        error("Failed to perform rectangular association analysis");
    }
    // write the results

    std::string assoc_outf = outf + assoc_suffix;
    htsFile* wf = hts_open(assoc_outf.c_str(), assoc_outf.substr(assoc_outf.length() - 3).compare(".gz") == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output file %s for writing", outf.c_str());
    }
    // write the header line
    hprintf(wf, "#CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tN_RR\tN_RA\tN_AA\tTEST");
    for(int32_t i = 0; i < input.pheno_matrix.pheno_ids.size(); ++i) {
        hprintf(wf, "\tBETA.%s\tSE.%s\tTSTAT.%s\tLOG10P.%s",
                input.pheno_matrix.pheno_ids[i].c_str(),
                input.pheno_matrix.pheno_ids[i].c_str(),
                input.pheno_matrix.pheno_ids[i].c_str(),
                input.pheno_matrix.pheno_ids[i].c_str());
    }
    hprintf(wf, "\n");
    for(int32_t j=0; j < n_loaded_vars; ++j) {
        const cpra_t& cpra = input.geno_chunk.v_cpra[j];
        const var_cnt_t& vcnt = input.geno_chunk.var_cnts[j];
        double a1freq = (double)(vcnt.ac) / (double)(vcnt.an);
        hprintf(wf, "%s\t%d\t%s\t%s\t%s\t%.6g\t%d\t%d\t%d\t%d\tLinear",
            cpra.chrom.c_str(),
            cpra.pos,
            cpra.to_string().c_str(),
            cpra.ref.c_str(),
            cpra.alts.c_str(),
            a1freq,
            vcnt.gcs[0] + vcnt.gcs[1] + vcnt.gcs[2],
            vcnt.gcs[0],
            vcnt.gcs[1],
            vcnt.gcs[2]);
        const std::vector<slr_sumstat_t>& var_results = rect_results[j];
        for(int32_t k=0; k < (int32_t)input.pheno_matrix.pheno_ids.size(); ++k) {
            const slr_sumstat_t& ss = var_results[k];
            hprintf(wf, "\t%.6g\t%.6g\t%.6g\t%.6g",
                ss.beta,
                ss.se,
                ss.tstat,
                ss.log10p);
        }
        hprintf(wf, "\n");
    }
    hts_close(wf);

    // ---- SuSiE fine-mapping (optional) --------------------------------------
    // Fine-mapping reuses the same genotype chunk and covariate-adjusted phenotype
    // matrix as the marginal association above, so the per-variant marginal BETA /
    // LOG10P (from rect_results) are matched directly into the credible-set output.
    if ( run_susie ) {
        const int32_t n_pheno = (int32_t)input.pheno_matrix.pheno_ids.size();
        const int32_t n_vars = input.geno_chunk.n_variants;
        if ( n_vars == 0 ) {
            notice("Skipping SuSiE: no variants loaded in the region");
        }
        else {
            susie::SusieOptions sopt;
            sopt.L = susie_L;
            sopt.max_iter = susie_max_iter;
            sopt.tol = susie_tol;
            sopt.standardize = !susie_no_standardize;

            // Unmappable-effects model. "inf"/"ash" force PIP-based convergence
            // (they have no well-defined ELBO). Match susieR's default inf/ash
            // tolerance of 1e-4 when the user left --susie-tol at its default.
            bool run_inf = false;
            bool run_ash = false;
            if ( unmappable_effects == "inf" ) {
                sopt.unmappable_effects = susie::SusieOptions::INF;
                sopt.convergence_method = susie::SusieOptions::PIP;
                if ( susie_tol == 1e-3 ) sopt.tol = 1e-4;
                run_inf = true;
            }
            else if ( unmappable_effects == "ash" ) {
                sopt.unmappable_effects = susie::SusieOptions::ASH;
                sopt.convergence_method = susie::SusieOptions::PIP;
                if ( susie_tol == 1e-3 ) sopt.tol = 1e-4;
                run_ash = true;
                // Reduction-test hooks: --ash-fix-pi / --ash-fix-sa2.
                auto parse_csv = [](const std::string& s) {
                    std::vector<double> out;
                    size_t i = 0;
                    while (i < s.size()) {
                        size_t j = s.find(',', i);
                        std::string tok = s.substr(i, j == std::string::npos ? std::string::npos : j - i);
                        if (!tok.empty()) out.push_back(atof(tok.c_str()));
                        if (j == std::string::npos) break;
                        i = j + 1;
                    }
                    return out;
                };
                if ( !ash_fix_pi_str.empty() ) sopt.ash_fix_pi = parse_csv(ash_fix_pi_str);
                if ( !ash_fix_sa2_str.empty() ) sopt.ash_fix_sa2 = parse_csv(ash_fix_sa2_str);
                if ( !sopt.ash_fix_pi.empty() && !sopt.ash_fix_sa2.empty() &&
                     sopt.ash_fix_pi.size() != sopt.ash_fix_sa2.size() ) {
                    error("--ash-fix-pi (K=%d) and --ash-fix-sa2 (K=%d) must have the same length",
                          (int)sopt.ash_fix_pi.size(), (int)sopt.ash_fix_sa2.size());
                }
                if ( !sopt.ash_fix_pi.empty() ) sopt.ash_K = (int)sopt.ash_fix_pi.size();
            }
            else if ( unmappable_effects != "none" ) {
                error("--unmappable-effects must be 'none', 'inf', or 'ash' (got '%s')", unmappable_effects.c_str());
            }
            const bool has_theta = run_inf || run_ash;

            // region label shared across all rows (CHR_BEG_END)
            char region_buf[256];
            snprintf(region_buf, sizeof(region_buf), "%s_%d_%d",
                     region_cbe.chrom.c_str(), region_cbe.beg1, region_cbe.end0);
            const std::string region_str(region_buf);

            // credible-set output: one row per variant that belongs to a credible set
            std::string cs_path = outf + susie_cs_suffix;
            htsFile* wcs = hts_open(cs_path.c_str(), cs_path.substr(cs_path.length() - 3).compare(".gz") == 0 ? "wz" : "w");
            if ( wcs == NULL )
                error("Cannot open SuSiE credible-set output file %s", cs_path.c_str());
            hprintf(wcs, "#trait\tvariant\tpip\tcs_id\talpha\tregion\tn_region_vars\tcs_size\tcs_lbf\tvar_lbf\tmu\tmu2\taf\tn\tbeta\tse\tlog10p\n");

            // optional per-variant LBF output: one row per variant, one column per single effect
            const int32_t L_eff = std::max(1, std::min(susie_L, n_vars));
            htsFile* wlbf = NULL;
            std::string lbf_path;
            if ( output_lbf ) {
                lbf_path = outf + susie_lbf_suffix;
                wlbf = hts_open(lbf_path.c_str(), lbf_path.substr(lbf_path.length() - 3).compare(".gz") == 0 ? "wz" : "w");
                if ( wlbf == NULL )
                    error("Cannot open SuSiE LBF output file %s", lbf_path.c_str());
                hprintf(wlbf, "#trait\tvariant\tregion\tn_region_vars\taf\tpip");
                // inf/ash: report per-variant unmappable-effect posterior mean
                if ( has_theta ) hprintf(wlbf, "\ttheta");
                for(int32_t l = 0; l < L_eff; ++l) hprintf(wlbf, "\tlbf.L%d", l + 1);
                hprintf(wlbf, "\n");
            }

            notice("Running SuSiE fine-mapping for %d phenotype(s) over %d variants", n_pheno, n_vars);
            for(int32_t k = 0; k < n_pheno; ++k) {
                const std::string& trait = input.pheno_matrix.pheno_ids[k];
                Eigen::VectorXd y = input.pheno_matrix.pheno_mat.col(k);
                susie::SusieResult res = susie::simple_susie_without_missing(
                    y, input.geno_chunk.geno_mat, sopt, susie_coverage, susie_min_abs_corr);

                const char* tag = run_ash ? "SuSiE-ash" : (run_inf ? "SuSiE-inf" : "SuSiE");
                if ( has_theta )
                    notice("%s trait %s: %d credible set(s), niter=%d, converged=%d, sigma2=%.4g, tau2=%.4g",
                           tag, trait.c_str(), (int)res.cs.size(), res.fit.niter,
                           (int)res.fit.converged, res.fit.sigma2, res.fit.tau2);
                else
                    notice("%s trait %s: %d credible set(s), niter=%d, converged=%d, sigma2=%.4g",
                           tag, trait.c_str(), (int)res.cs.size(), res.fit.niter,
                           (int)res.fit.converged, res.fit.sigma2);

                // one row per variant in each credible set
                for(int32_t c = 0; c < (int32_t)res.cs.size(); ++c) {
                    const susie::CredibleSet& cs = res.cs[c];
                    const int32_t l = cs.effect_index;      // single effect defining this CS
                    const int32_t cs_id = c + 1;            // 1-based CS id within the trait
                    const int32_t cs_size = (int32_t)cs.variables.size();
                    for(int32_t m = 0; m < cs_size; ++m) {
                        const int32_t j = cs.variables[m];
                        const cpra_t& cpra = input.geno_chunk.v_cpra[j];
                        const var_cnt_t& vcnt = input.geno_chunk.var_cnts[j];
                        const slr_sumstat_t& ss = rect_results[j][k]; // marginal association
                        const double af = (double)vcnt.ac / (double)vcnt.an;
                        // per-variant log BF (natural log) under the single effect defining this CS
                        const double var_lbf = (l < (int32_t)res.fit.lbf_variable.rows())
                            ? res.fit.lbf_variable(l, j) : 0.0;
                        hprintf(wcs, "%s\t%s\t%.6g\t%d\t%.6g\t%s\t%d\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\t%d\t%.6g\t%.6g\t%.6g\n",
                            trait.c_str(),
                            cpra.to_string().c_str(),   // variant C:P:R:A
                            res.fit.pip(j),             // marginal PIP
                            cs_id,
                            res.fit.alpha(l, j),        // posterior prob within this CS
                            region_str.c_str(),
                            n_vars,                     // number of variants tested in the region
                            cs_size,
                            cs.lbf,                     // log BF of the single effect defining this CS
                            var_lbf,                    // log BF of this variant under that effect
                            res.fit.mu(l, j),           // posterior mean | inclusion
                            res.fit.mu2(l, j),          // posterior 2nd moment | inclusion
                            af,
                            ss.n_obs,                   // sample size
                            ss.beta,                    // marginal BETA
                            ss.se,                      // marginal SE
                            ss.log10p);                 // marginal LOG10P
                    }
                }

                // per-variant LBF for every variant (optional)
                if ( output_lbf ) {
                    const int32_t Lrows = (int32_t)res.fit.lbf_variable.rows();
                    for(int32_t j = 0; j < n_vars; ++j) {
                        const cpra_t& cpra = input.geno_chunk.v_cpra[j];
                        const var_cnt_t& vcnt = input.geno_chunk.var_cnts[j];
                        const double af = (double)vcnt.ac / (double)vcnt.an;
                        // marginal PIP = 1 - prod_l(1 - alpha_lj), reported for every
                        // variant regardless of credible-set membership
                        hprintf(wlbf, "%s\t%s\t%s\t%d\t%.6g\t%.6g", trait.c_str(), cpra.to_string().c_str(), region_str.c_str(), n_vars, af, res.fit.pip(j));
                        // inf/ash: unmappable-effect posterior mean (standardized-X scale, like susieR fit$theta)
                        if ( has_theta ) hprintf(wlbf, "\t%.6g", (j < (int32_t)res.fit.theta.size()) ? res.fit.theta(j) : 0.0);
                        for(int32_t l = 0; l < L_eff; ++l) {
                            const double v = (l < Lrows) ? res.fit.lbf_variable(l, j) : 0.0;
                            hprintf(wlbf, "\t%.6g", v);
                        }
                        hprintf(wlbf, "\n");
                    }
                }
            }
            hts_close(wcs);
            if ( wlbf ) hts_close(wlbf);
            notice("SuSiE credible sets written to %s", cs_path.c_str());
            if ( output_lbf ) notice("SuSiE per-variant LBF written to %s", lbf_path.c_str());
        }
    }

    notice("Analysis finished. Total variants processed: %lld", (long long)input.geno_chunk.n_variants);
    return 0;
}