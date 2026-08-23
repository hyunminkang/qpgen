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

int32_t cmd_rect_assoc(int32_t argc, char **argv)
{
    std::string pgenf;      // PLINK2 genotype file
    std::string psamf;      // PLINK2 sample file
    std::string pivarf;     // PLINK2 indexed pvar file, compatible with PgenIdxReader
    std::string pgenlistf;  // CHROM BEG END PGEN PSAM PIVAR containing the list of region-specific BED files
    std::string phef;
    std::string covf;
    std::string varlistf;   // List of variants to be tested
    std::string phelistf;   // List of phenotypes to be tested
    std::string outf;
    std::string samplef;
    std::string pheno_format("regenie");
    std::string cov_format("regenie");
    int32_t jump_thres_bp = 1000000;
    int32_t max_chunk_vars = 100; // Maximum number of variants to store at once in memory
    int32_t icol_pivar_idx = 9;
    int32_t icol_pheno_id = 1;
    int32_t icol_cov_id = 1;
    double min_maf = 1e-10;
    double min_mac = 1.0;
    bool rint_before_adj = false; // Perform rank-based inverse normal transformation before covariate adjustment
    bool rint_after_adj = false; // Perform rank-based inverse normal transformation after covariate adjustment

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("pgen", &pgenf, "Input PLINK2 genotype file")
    LONG_STRING_PARAM("psam", &psamf, "Input PLINK2 sample file")
    LONG_STRING_PARAM("pivar", &pivarf, "Input PLINK2 index pvar file (bgzipped and tabix)")
    LONG_STRING_PARAM("pgen-list", &pgenlistf, "Input file containing the list of region-specific BED files")
    LONG_STRING_PARAM("pheno", &phef, "Input phenotype matrix")
    LONG_STRING_PARAM("sample", &samplef, "Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files")
    LONG_STRING_PARAM("cov", &covf, "Input covariate matrix (optional)")
    LONG_STRING_PARAM("var-list", &varlistf, "Input file containing the list of variants to be tested")
    LONG_STRING_PARAM("pheno-list", &phelistf, "Input file containing the list of phenotypes to be tested. If not provided, all phenotypes in the phenotype matrix will be tested.")
    LONG_STRING_PARAM("pheno-format", &pheno_format, "Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("cov-format", &cov_format, "Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")

    LONG_PARAM_GROUP("Analysis options", NULL)
    LONG_PARAM("rint-before-adj", &rint_before_adj, "Perform rank-based inverse normal transformation before covariate adjustment (default: false)")
    LONG_PARAM("rint-after-adj", &rint_after_adj, "Perform rank-based inverse normal transformation after covariate adjustment (default: false)")
    LONG_INT_PARAM("max-chunk-vars", &max_chunk_vars, "Maximum number of variants to store at once in memory (default: 100)")
    LONG_DOUBLE_PARAM("min-maf", &min_maf, "Minimum minor allele frequency to include a variant in the analysis (default: 0.0)")
    LONG_DOUBLE_PARAM("min-mac", &min_mac, "Minimum minor allele count to include a variant in the analysis (default: 0.0)")

    LONG_PARAM_GROUP("Auxiliary options", NULL)
    LONG_INT_PARAM("jump-thres-bp", &jump_thres_bp, "Jump threshold in base pairs for the variant index (default: 1000000)")
    LONG_INT_PARAM("icol-pivar-idx", &icol_pivar_idx, "1-based column index for the variant ID in the pvar file (default: 9)")
    LONG_INT_PARAM("icol-pheno-id", &icol_pheno_id, "1-based column index for the phenotype ID in the phenotype file (default: 1)")
    LONG_INT_PARAM("icol-cov-id", &icol_cov_id, "1-based column index for the covariate ID in the phenotype file (default: 1)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // check required arguments
    if ( phef.empty() ) {
        error("Phenotype file (--pheno) is required");
    }
    if ( varlistf.empty() ) {
        error("Variant list file (--var-list) is required");
    }
    if ( outf.empty() ) {
        error("Output prefix (--out) is required");
    }

    // set up the association input object and configure the parameters. This is the
    // same plumbing cmd_region_assoc uses, and gives us the covariate-adjusted phenotype
    // matrix plus the overlapping sample IDs we need to load genotypes against.
    ind_assoc_input input;
    input.set_jump_thres_bp(jump_thres_bp);
    input.set_icol_pivar_idx(icol_pivar_idx - 1); // convert to 0-based index
    input.set_rint_before_adj(rint_before_adj);
    input.set_rint_after_adj(rint_after_adj);
    if ( !samplef.empty() ) {
        input.set_subset_sample_file(samplef.c_str());
    }

    if ( !pgenlistf.empty() ) {
        if ( !pgenf.empty() || !pivarf.empty() || !psamf.empty() ) {
            error("When --pgen-list is provided, --pgen, --pivar, and --psam should not be provided");
        }
        input.process_pgenlist(pgenlistf.c_str(), phef.c_str(), pheno_format.c_str(), covf.c_str(), cov_format.c_str());
    }
    else if ( !pgenf.empty() && !pivarf.empty() && !psamf.empty() ) {
        input.process_single_pgen(pgenf.c_str(), pivarf.c_str(), psamf.c_str(), phef.c_str(), pheno_format.c_str(), covf.c_str(), cov_format.c_str());
    }
    else {
        error("Either --pgen-list or (--pgen, --pivar, and --psam) must be provided");
    }

    const int32_t n_overlapping_samples = (int32_t)input.overlapping_sample_ids.size();

    // load_pheno_cov_matrices already applies rint-before-adj / covariate adjustment /
    // rint-after-adj to input.pheno_matrix, so from here we just consume input.pheno_matrix.pheno_mat.
    PhenoMatrix& pheno_matrix = input.pheno_matrix;

    // load the trait list -- keep the permissive overlap semantics the original had:
    // phenotypes named in --pheno-list that are not in the pheno file are silently dropped
    if ( !phelistf.empty() ) {
        tsv_reader tr_phelist(phelistf.c_str());
        std::vector<std::string> phelist_ids;
        while( tr_phelist.read_line() ) {
            if ( tr_phelist.nfields != 1 ) {
                error("Invalid format phenotype list file %s in line %zu starting with %s. Must contain only 1 field", phelistf.c_str(), (int32_t)phelist_ids.size() + 1, tr_phelist.str_field_at(0) );
            }
            phelist_ids.push_back( tr_phelist.str_field_at(0) );
        }
        std::vector<std::string> overlapping_phe_ids;
        notice("Loaded %d phenotype IDs from %s to subset the phenotype matrix", (int32_t)phelist_ids.size(), phelistf.c_str());
        identify_overlapping_ids(phelist_ids, pheno_matrix.pheno_ids, overlapping_phe_ids);
        notice("%d overlapping phenotypes found between %d phenotypes in the phenotype matrix and %d phenotypes in %s", (int32_t)overlapping_phe_ids.size(), (int32_t)pheno_matrix.pheno_ids.size(), (int32_t)phelist_ids.size(), phelistf.c_str());
        if ( pheno_matrix.subset_pheno_ids(overlapping_phe_ids) != (int32_t)overlapping_phe_ids.size() ) {
            error("Failed to subset the phenotype matrix based on the provided phenotype list in %s", phelistf.c_str());
        }
    }

    // load the variant list
    tsv_reader tr_varlist(varlistf.c_str()); // varlist can be either [CHROM]:[POS]:[REF]:[ALT] format or tab-delimited
    std::vector<std::string> varlist_cpra;
    notice("Reading the variant list from %s", varlistf.c_str());
    while( tr_varlist.read_line() ) {
        if ( tr_varlist.nfields == 1 ) {
            varlist_cpra.push_back(tr_varlist.str_field_at(0));
        }
        else if ( tr_varlist.nfields >= 4 ) {
            char buf[65536];
            snprintf(buf, sizeof(buf), "%s:%s:%s:%s",
                tr_varlist.str_field_at(0),
                tr_varlist.str_field_at(1),
                tr_varlist.str_field_at(2),
                tr_varlist.str_field_at(3));
            varlist_cpra.push_back(buf);
        }
        else {
            error("Invalid format in the variant list file %s at line %zu starting with %s", varlistf.c_str(), (int32_t)varlist_cpra.size() + 1, tr_varlist.str_field_at(0));
        }
    }

    // whether covariate residualization of the genotype chunks is required (Frisch-Waugh-Lovell).
    // load_pheno_cov_matrices adjusted the phenotype matrix; to estimate partial effects we must
    // adjust the genotype matrix against the same covariates as well.
    const bool have_cov = !covf.empty() && input.cov_matrix.pheno_mat.cols() > 0;

    // open the output file gz or plain based on the extension
    htsFile* wf = hts_open(outf.c_str(), outf.substr(outf.length() - 3).compare(".gz") == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output file %s for writing", outf.c_str());
    }
    // write the header line
    hprintf(wf, "#CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tN_RR\tN_RA\tN_AA\tTEST");
    for(int32_t i = 0; i < pheno_matrix.pheno_ids.size(); ++i) {
        hprintf(wf, "\tBETA.%s\tSE.%s\tTSTAT.%s\tLOG10P.%s",
                pheno_matrix.pheno_ids[i].c_str(),
                pheno_matrix.pheno_ids[i].c_str(),
                pheno_matrix.pheno_ids[i].c_str(),
                pheno_matrix.pheno_ids[i].c_str());
    }
    hprintf(wf, "\n");

    // process each chunk
    for(int32_t i=0; i < varlist_cpra.size(); i += max_chunk_vars) {
        bool geno_has_missing = false;
        int32_t chunk_size = (int32_t)varlist_cpra.size() - i > max_chunk_vars ? max_chunk_vars : (int32_t)varlist_cpra.size() - i;
        Eigen::MatrixXd geno_mat(n_overlapping_samples, chunk_size);
        Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> geno_mask(n_overlapping_samples, chunk_size);
        int32_t n_skipped = 0;
        notice("Loading genotype data for chunk %d with %d variants", i / max_chunk_vars + 1, chunk_size);
        std::vector<std::string> chunk_cpras;
        std::vector<var_cnt_t> chunk_var_cnts;
        for(int32_t j=0; j < chunk_size; ++j) {
            int32_t gcs[3] = {0, 0, 0};
            if ( !input.mpr.read_pivar(varlist_cpra[i + j].c_str()) ) {
                notice("Skipping variant %s, which is not found in the genotype data", varlist_cpra[i + j].c_str());
                ++n_skipped;
                continue;
            }
            else {
                // load the genotype data
                if ( !input.mpr.get_genos() ) {
                    notice("Cannot load genotype data for variant %s. Skipping", varlist_cpra[i + j].c_str());
                    ++n_skipped;
                    continue;
                }
                int32_t jv = j - n_skipped;
                const std::vector<int32_t>& int_buf = input.mpr.get_int_buf();
                const double* dbl_buf = input.mpr.get_dbl_buf();
                int32_t an = 0;
                double ac = 0;
                double info = 0;
                if ( input.mpr.is_dosage_present() ) {
                    if ( dbl_buf == NULL) {
                        dbl_buf = input.mpr.get_dbl_buf();
                    }
                    double sumsq = 0;
                    for(int32_t k =0; k < n_overlapping_samples; ++k) {
                        double ds = 2.0 - dbl_buf[k];
                        geno_mat(k, jv) = ds;
                        geno_mask(k, jv) = true; // not missing
                        ac += ds;
                        an += 2;
                        ++gcs[(ds < 0.5) ? 0 : ( (ds < 1.5) ? 1 : 2 )];
                        sumsq += (ds * ds);
                    }
                    double af = ac / (double)an;
                    for(int32_t k =0; k < n_overlapping_samples; ++k) {
                        geno_mat(k, jv) -= (2 * af); // center to zero
                    }

                    if ( ac < min_mac || af < min_maf || (1.0 - af) < min_maf  || an - ac < min_mac ) {
                        // skip rare variants below the MAF or MAC threshold
                        ++n_skipped;
                        continue;
                    }

                    // E(Var(g)) = 2 * af * (1-af) = 2 * ac * ( an - ac ) / an / an
                    // Var(g) = sumsq / (an / 2) - 4 * ac * ac / an / an
                    // ratio = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) )
                    info = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) );
                }
                else {
                    for(int32_t k = 0; k < n_overlapping_samples; ++k) {
                        switch(int_buf[k]) { // make sure to convert 1-based index to 0-based
                        case 0:
                            an += 2;
                            ac += 2;
                            ++gcs[2];
                            break;
                        case 1:
                            an += 2;
                            ++ac;
                            ++gcs[1];
                            break;
                        case 2:
                            an += 2;
                            ++gcs[0];
                            break;
                        }
                    }

                    double af = ac / (double)an;
                    if ( ac < min_mac || af < min_maf || (1.0 - af) < min_maf  || an - ac < min_mac ) {
                        // skip rare variants below the MAF or MAC threshold
                        ++n_skipped;
                        continue;
                    }

                    double mean = (double)ac / (double)an * 2.0;
                    // Exp(Var(g)) = af * (1-af) * 2 = mean * (2 - mean) / 2;
                    // Var(g) = EX^2 - EX^2 = (4 * n_2 + 1 * n_1)/n - mean^2
                    info = ((4.0 * gcs[2] + gcs[1])/(an/2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
                    for(int32_t k = 0; k < n_overlapping_samples; ++k) {
                        switch(int_buf[k]) {
                        case 0:
                            geno_mat(k, jv) = 2.0 - mean; // homalt
                            geno_mask(k, jv) = true; // not missing
                            break;
                        case 1:
                            geno_mat(k, jv) = 1.0 - mean; // het
                            geno_mask(k, jv) = true; // not missing
                            break;
                        case 2:
                            geno_mat(k, jv) = 0.0 - mean; // homref
                            geno_mask(k, jv) = true; // not missing
                            break;
                        default:
                            geno_mat(k, jv) = 0; // missing - mean imputation
                            geno_mask(k, jv) = false; // missing
                            geno_has_missing = true;
                            break;
                        }
                    }
                }
                chunk_cpras.push_back(varlist_cpra[i + j]);
                chunk_var_cnts.push_back( var_cnt_t( an, ac, gcs[0], gcs[1], gcs[2] ) );
            }
        }
        int32_t new_chunk_size = chunk_size - n_skipped;
        geno_mat.conservativeResize(Eigen::NoChange, new_chunk_size);
        geno_mask.conservativeResize(Eigen::NoChange, new_chunk_size);

        // Residualize the genotypes against the same covariates used to adjust the phenotype
        // (Frisch-Waugh-Lovell). See load_genotype_chunk() in assoc_utils.cpp -- adjusting only
        // the phenotype leaves genotype variance collinear with the covariates in the design,
        // which shrinks beta and the test statistic. Applied here to match cmd_region_assoc.
        if ( have_cov && new_chunk_size > 0 ) {
            if ( input.cov_matrix.pheno_mat.rows() != n_overlapping_samples ) {
                error("Covariate matrix has %d rows but %d overlapping samples were expected",
                      (int32_t)input.cov_matrix.pheno_mat.rows(), n_overlapping_samples);
            }
            geno_mat = pheno_adj_cov_nxt_without_missing(geno_mat, input.cov_matrix.pheno_mat);
        }

        // perform rectangular association analysis
        std::vector<std::vector<slr_sumstat_t> > rect_results;
        notice("Performing rectangular association analysis for %d phenotypes and %d variants after skipping %d variants", (int32_t)pheno_matrix.pheno_ids.size(), new_chunk_size, n_skipped);
        if ( !simple_rect_regression_without_missing(
                pheno_matrix.pheno_mat,
                geno_mat,
                rect_results) ) {
            error("Failed to perform rectangular association analysis for chunk %d", i / max_chunk_vars + 1);
        }

        // write the results
        for(int32_t j=0; j < new_chunk_size; ++j) {
            cpra_t cpra(chunk_cpras[j].c_str());
            const var_cnt_t& vcnt = chunk_var_cnts[j];
            double a1freq = (double)(vcnt.ac) / (double)(vcnt.an);
            hprintf(wf, "%s\t%d\t%s\t%s\t%s\t%.6g\t%d\t%d\t%d\t%d\tLinear",
                cpra.chrom.c_str(),
                cpra.pos,
                chunk_cpras[j].c_str(),
                cpra.ref.c_str(),
                cpra.alts.c_str(),
                a1freq,
                vcnt.an,
                vcnt.gcs[0],
                vcnt.gcs[1],
                vcnt.gcs[2]);
            const std::vector<slr_sumstat_t>& var_results = rect_results[j];
            for(int32_t k=0; k < (int32_t)pheno_matrix.pheno_ids.size(); ++k) {
                const slr_sumstat_t& ss = var_results[k];
                hprintf(wf, "\t%.6g\t%.6g\t%.6g\t%.6g",
                    ss.beta,
                    ss.se,
                    ss.tstat,
                    ss.log10p);
            }
            hprintf(wf, "\n");
        }
    }
    hts_close(wf); // close the output file
    notice("Analysis finished");
    return 0;
}
