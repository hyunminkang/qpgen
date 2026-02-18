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

int32_t cmd_pair_assoc(int32_t argc, char **argv)
{
    std::string pgenf;      // PLINK2 genotype file
    std::string psamf;      // PLINK2 sample file
    std::string pivarf;     // PLINK2 indexed pvar file, compatible with PgenIdxReader
    std::string pgenlistf;      // CHROM BEG END PGEN PSAM PIVAR containing the list of region-specific BED files
    std::string phef;
    std::string covf;
    std::string pairf;
    std::string outf;
    std::string samplef;
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
    std::string colname_pair_trait("trait"); // column name for the trait ID in the pair file
    std::string colname_pair_variant("variant"); // column name of the variant ID in the pair file
    std::string colname_pair_region("region"); // column name of the region string in the pair file
    bool skip_rint = false; // skip tests based on rank-based inverse normal transformation

    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("pgen", &pgenf, "Input PLINK2 genotype file")
    LONG_STRING_PARAM("psam", &psamf, "Input PLINK2 sample file")
    LONG_STRING_PARAM("pivar", &pivarf, "Input PLINK2 index pvar file (bgzipped and tabix)")
    LONG_STRING_PARAM("pgen-list", &pgenlistf, "Input file containing CHROM BEG END PGEN PSAM PIVAR")
    LONG_STRING_PARAM("pheno", &phef, "Input phenotype matrix")
    LONG_STRING_PARAM("pairs", &pairf, "Input file containing [trait_id] [variant_id] pairs to be tested")
    LONG_STRING_PARAM("sample", &samplef, "Input file containing sample IDs to be used. Useful when different IDs are used in pgen and pheno files")
    LONG_STRING_PARAM("cov", &covf, "Input covariate matrix (optional)")
    LONG_STRING_PARAM("pheno-format", &pheno_format, "Format of the phenotype file (default: 'regenie'). Options: 'regenie', 'tensorqtl', 'tsv-sample-col', 'tsv-sample-row'")
    LONG_STRING_PARAM("cov-format", &cov_format, "Format of the covariate file (default: 'regenie'). Options: 'regenie', 'tsv-sample-col', 'tsv-sample-row'")


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
    LONG_STRING_PARAM("colname-pair-trait", &colname_pair_trait, "Column name of the trait ID in the pait file")
    LONG_STRING_PARAM("colname-pair-variant", &colname_pair_variant, "Column name of the variant ID in the pair file")
    LONG_STRING_PARAM("colname-pair-region", &colname_pair_region, "Column name of the region string in the pair file")
    LONG_PARAM("skip-rint", &skip_rint, "Skip tests based on rank-based inverse normal transformation (default: false)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // open the Plink file
    MultiPgenIdxReader mpr;
    if ( !pgenlistf.empty() ) { // list is provided
        //notice("foo");
        if ( pgenf.empty() && pivarf.empty() && psamf.empty() ) {
            if ( !mpr.prep_pgen_list(pgenlistf.c_str()) ) {
                error("Failed to prepare pgen files with the following list file: %s", pgenlistf.c_str());
            }
        }
        else {
            error("When --pgen-list is provided, --pgen, --pivar, and --psam should not be provided");
        }
    }
    else {
        //notice("bar");
        if ( pgenf.empty() || pivarf.empty() || psamf.empty() ) {
            error("When --pgen-list is not provided, --pgen, --pivar, and --psam should be provided");
        }
        if ( !mpr.set_single_chunk_pgen(pgenf.c_str(), pivarf.c_str(), psamf.c_str()) ) {
            error("Failed to add pgen files with the following files:\n%s\n%s\n%s", pgenf.c_str(), pivarf.c_str(), psamf.c_str());
        }
    }
    mpr.set_jump_thres_bp(jump_thres_bp);
    mpr.set_icol_pivar_idx(icol_pivar_idx - 1); // convert to 0-based index

    notice("Prepared to load PLINK2 genotype data with %d samples", mpr.get_all_sample_count());

    // load the sample IDs
    std::vector<std::string> subset_sample_ids;
    if ( !samplef.empty() ) {
        notice("Loading sample IDs from %s", samplef.c_str());
        tsv_reader tr_sample(samplef.c_str());
        while(tr_sample.read_line()) {
            subset_sample_ids.push_back(tr_sample.str_field_at(0));
        }
        notice("Loaded %d sample IDs from %s to subset the genotype and phenotype data", (int32_t)subset_sample_ids.size(), samplef.c_str());
    }

    // find overlapping sample IDs
    const std::vector<plink_samp_t>& geno_all_samps = mpr.get_all_samples();
    std::vector<std::string> geno_all_samp_ids;
    for(int32_t i=0; i < geno_all_samps.size(); ++i) {
        geno_all_samp_ids.push_back(geno_all_samps[i].indID);
    }

    notice("Loading phenotype matrix from %s", phef.c_str());

    // load the phenotype matrix
    PhenoMatrix pheno_matrix;
    if ( !pheno_matrix.load_pheno_matrix(phef.c_str(), pheno_format.c_str()) ) {
        error("Failed to load the phenotype matrix from file %s", phef.c_str());
    }

    notice("Loaded phenotype matrix with %d samples and %d phenotypes from %s", (int32_t)pheno_matrix.samp_ids.size(), (int32_t)pheno_matrix.pheno_ids.size(), phef.c_str());

    PhenoMatrix cov_matrix;
    if ( !covf.empty() ) {
        notice("Loading covariate matrix from %s", covf.c_str());
        if ( !cov_matrix.load_pheno_matrix(covf.c_str(), cov_format.c_str()) ) {
            error("Failed to load the covariate matrix from file %s", covf.c_str());
        }
        notice("Loaded covariate matrix with %d samples and %d covariates from %s", (int32_t)cov_matrix.samp_ids.size(), (int32_t)cov_matrix.pheno_ids.size(), covf.c_str());
    }

    // load the covariate matrix
    std::vector<std::string> overlapping_sample_ids;
    if ( !subset_sample_ids.empty() ) {
        std::vector<std::string> temp1_ids;
        std::vector<std::string> temp2_ids;
        identify_overlapping_ids(subset_sample_ids, geno_all_samp_ids, temp1_ids);
        if ( covf.empty() ) {
            identify_overlapping_ids(temp1_ids, pheno_matrix.samp_ids, overlapping_sample_ids);
        }
        else {
            identify_overlapping_ids(temp1_ids, cov_matrix.samp_ids, temp2_ids);
            identify_overlapping_ids(temp2_ids, pheno_matrix.samp_ids, overlapping_sample_ids);
        }
    }
    else {
        if ( covf.empty() ) {
            identify_overlapping_ids(geno_all_samp_ids, pheno_matrix.samp_ids, overlapping_sample_ids);
        }
        else {
            std::vector<std::string> temp1_ids;
            identify_overlapping_ids(geno_all_samp_ids, cov_matrix.samp_ids, temp1_ids);
            identify_overlapping_ids(temp1_ids, pheno_matrix.samp_ids, overlapping_sample_ids);
        }
    }
    notice("%zu overlapping samples found among sample, genotype, phenotype, and covariate files", (int32_t)overlapping_sample_ids.size());

    if ( pheno_matrix.samp_ids.size() != overlapping_sample_ids.size() ) {
        notice("Subsetting the phenotype matrix to the overlapping samples");
        pheno_matrix.subset_sample_ids(overlapping_sample_ids);
    }
    if ( !covf.empty() ) {
        if ( cov_matrix.samp_ids.size() != overlapping_sample_ids.size() ) {
            notice("Subsetting the covariate matrix to the overlapping samples");
            cov_matrix.subset_sample_ids(overlapping_sample_ids);
        }
    }
    if ( mpr.get_all_sample_count() != overlapping_sample_ids.size() ) {
        notice("Subsetting the genotype data to %zu overlapping samples", (int32_t)overlapping_sample_ids.size());
        mpr.subset_sample_ids(overlapping_sample_ids);
    }
    int32_t n_overlapping_samples = (int32_t)overlapping_sample_ids.size();

    // adjust phenotype matrix by covariates
    // if ( rint_before_adj ) {
    //     notice("Performing rank-based inverse normal transformation for all phenotypes before covariate adjustment");
    //     if ( !pheno_matrix.has_missing ) {
    //         pheno_matrix.pheno_mat = rint_matrix_without_missing(pheno_matrix.pheno_mat);
    //     }
    //     else {
    //         error("Rank-based inverse normal transformation adjustment is currently only supported for phenotype matrices with missing values");
    //     }
    // }

    // perform covariate adjustment
    if ( !covf.empty() ) {
        notice("Adjusting phenotypes by covariates using linear regression");
        if ( pheno_matrix.has_missing || cov_matrix.has_missing ) {
            error("Covariate adjustment is currently only supported for phenotype and covariate matrices without missing values");
        }
        else {
            pheno_matrix.pheno_mat = pheno_adj_cov_nxt_without_missing(pheno_matrix.pheno_mat, cov_matrix.pheno_mat);
        }
    }

    // if ( rint_after_adj ) {
    //     notice("Performing rank-based inverse normal transformation for all phenotypes after covariate adjustment");
    //     if ( !pheno_matrix.has_missing ) {
    //         pheno_matrix.pheno_mat = rint_matrix_without_missing(pheno_matrix.pheno_mat);
    //     }
    //     else {
    //         error("Rank-based inverse normal transformation adjustment is currently only supported for phenotype matrices without missing values");
    //     }
    // }

    // PhenoMatrix pheno_matrix_rint;
    // if ( !skip_rint ) { 
    //     pheno_matrix_rint = pheno_matrix; // copy constructor
    //     notice("Preparing rank-based inverse normal transformed phenotypes for all traits");
    //     if ( !pheno_matrix_rint.has_missing ) {
    //         pheno_matrix_rint.pheno_mat = rint_matrix_without_missing(pheno_matrix_rint.pheno_mat);
    //     }
    //     else {
    //         error("Rank-based inverse normal transformation adjustment is currently only supported for phenotype matrices without missing values");
    //     }
    // }

    // read the variant-trait pair lists
    tsv_reader tr_pair(pairf.c_str());
    std::map<int32_t, std::set<cpra_t> > pair_map;
    std::map<int32_t, std::set<cpra_t> >::iterator pair_map_it;
    std::map<int32_t, std::set<cbe_t> > region_map;
    std::map<int32_t, std::set<cbe_t> >::iterator region_map_it;

    notice("Reading the trait-variant/region pairs from %s", pairf.c_str());
    int32_t n_pairs = 0;

    // open the output file gz or plain based on the extension
    htsFile* wf = hts_open(outf.c_str(), outf.substr(outf.length() - 3).compare(".gz") == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output file %s for writing", outf.c_str());
    }
    // write the header line
    //hprintf(wf, "#TRAIT\tCHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tN_RR\tN_RA\tN_AA\tTEST");
    //hprintf(wf, "#TRAIT\tCHROM\tPOS\tID\tREF\tALT\tAF\tINFO\tN\tBETA\tSE\tTSTAT\tLOG10P");
    hprintf(wf, "#TRAIT\tCHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tN_RR\tN_RA\tN_AA\tINFO\tBETA\tSE\tTSTAT\tLOG10P");
    if ( !skip_rint ) {
        hprintf(wf, "\tBETA_RINT\tSE_RINT\tTSTAT_RINT\tLOG10P_RINT");
    }
    hprintf(wf, "\n");

    int32_t icol_pair_trait = -1;
    int32_t icol_pair_variant = -1;
    int32_t icol_pair_region = -1;
    bool region_mode = false;
    std::map<std::string, int32_t>::iterator phe_trait2idx_it;
    while( tr_pair.read_line() ) {
        if ( icol_pair_trait < 0 ) {
            for(int32_t i = 0; i < tr_pair.nfields; ++i) {
                if ( colname_pair_trait.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_trait = i;
                }
                else if ( colname_pair_variant.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_variant = i;
                    region_mode = false;
                }
                else if ( colname_pair_region.compare(tr_pair.str_field_at(i)) == 0 ) {
                    icol_pair_region = i;
                    region_mode = true;
                }
            }
            if ( icol_pair_trait < 0 ) {
                error("Cannot find column %s in the pair file %s", colname_pair_trait.c_str(), pairf.c_str());
            }
            if ( icol_pair_variant < 0 && icol_pair_region < 0 ) {
                error("Either %s or %s must exist in the pair file %s", colname_pair_variant.c_str(), colname_pair_region.c_str(), pairf.c_str());
            }
            if ( icol_pair_variant >= 0 && icol_pair_region >= 0 ) {
                error("Both %s and %s cannot exist in the pair file %s. Please provide either variant ID or region string", colname_pair_variant.c_str(), colname_pair_region.c_str(), pairf.c_str());
            }
        }
        else {
            if ( icol_pair_variant >= 0 ) {
                const char* var_id = tr_pair.str_field_at(icol_pair_variant);
                const char* phe_id = tr_pair.str_field_at(icol_pair_trait);
                phe_trait2idx_it = pheno_matrix.pheno_id2idx.find(phe_id);
                if ( phe_trait2idx_it == pheno_matrix.pheno_id2idx.end() ) {
                    notice("Skipping phenotype %s, which is not observed in %s", phe_id, phef.c_str());
                }
                else {
                    int32_t phe_idx = phe_trait2idx_it->second;
                    cpra_t cpra(var_id);
                    if ( pair_map[phe_idx].insert(cpra).second ) {
                        ++n_pairs;
                    }
                    else {
                        notice("Skipping duplicate pair %s\t%s", var_id, phe_id);
                    }
                }
            }
            else {
                const char* phe_id = tr_pair.str_field_at(icol_pair_trait);
                const char* region_str = tr_pair.str_field_at(icol_pair_region);
                phe_trait2idx_it = pheno_matrix.pheno_id2idx.find(phe_id);
                if ( phe_trait2idx_it == pheno_matrix.pheno_id2idx.end() ) {
                    notice("Skipping phenotype %s, which is not observed in %s", phe_id, phef.c_str());
                }
                else {
                    int32_t phe_idx = phe_trait2idx_it->second;
                    cbe_t cbe(region_str);
                    if ( region_map[phe_idx].insert(cbe).second ) {
                        ++n_pairs;
                    }
                    else {
                        notice("Skipping duplicate pair %s\t%s", region_str, phe_id);
                    }
                }
                //error("region-based pair association is not supported yet. Please provide the variant ID in the pair file");
            }
        }
    }
    notice("Finished reading %d trait-%s pairs across %zu traits", n_pairs, region_mode ? "region" : "variant", pair_map.size());

    if ( region_mode ) {
        for(region_map_it = region_map.begin(); region_map_it != region_map.end(); ++region_map_it) {
            int32_t phe_idx = region_map_it->first;
            const std::set<cbe_t>& regions = region_map_it->second;
            if ( regions.empty() ) continue;

            notice("Processing %zu regions for phenotype %s", regions.size(), pheno_matrix.pheno_ids[phe_idx].c_str());
            Eigen::VectorXd phe_vec(n_overlapping_samples);
            phe_vec = pheno_matrix.pheno_mat.col(phe_idx);
            //phe_vec = phe_mat_adj.row(phe_idx); // get the phenotype vector for the current trait

            // process one region at a time. The variant may overlap, but let's not worry about it for now
            for(std::set<cbe_t>::const_iterator regions_it = regions.begin(); regions_it != regions.end(); ++regions_it) {
                const cbe_t& region = *regions_it;
                notice("Processing region %s for phenotype %s", region.to_string().c_str(), pheno_matrix.pheno_ids[phe_idx].c_str());
                //Eigen::Vector<bool, Eigen::Dynamic> phe_mask_vec = phe_mask.row(phe_idx); // get the missing values for the current trait
                Eigen::Vector<bool, Eigen::Dynamic> phe_mask_vec = pheno_matrix.pheno_mask.col(phe_idx); // get the missing values for the current trait

                bool phe_has_missing = !phe_mask_vec.all();
                Eigen::VectorXd phe_rint_vec;
                if ( ! skip_rint ) {
                    // perform rank-based inverse normal transformation
                    notice("Performing rank-based inverse normal transformation for phenotype %s", pheno_matrix.pheno_ids[phe_idx].c_str());
                    if ( phe_has_missing ) {
                        phe_rint_vec = rint_with_missing(phe_vec, pheno_matrix.pheno_mask.col(phe_idx));
                    }
                    else {
                        phe_rint_vec = rint_without_missing(phe_vec);
                    }
                    if ( phe_rint_vec.size() != n_overlapping_samples ) {
                        error("Rank-based inverse normal transformation failed for phenotype %s", pheno_matrix.pheno_ids[phe_idx].c_str());
                    }
                }
                else {
                    notice("Skipping rank-based inverse normal transformation for phenotype %s", pheno_matrix.pheno_ids[phe_idx].c_str());
                }

                if ( mpr.read_pos(region.chrom.c_str(), region.beg1) ) { // variant exists, start reading the region
                    std::vector<cpra_t> v_cpra;
                    std::vector<int32_t> ans;
                    std::vector<double> acs;
                    std::vector<int32_t> gc0s, gc1s, gc2s; 
                    std::vector<double> infos;
                    const std::vector<int32_t>& int_buf = mpr.get_int_buf();
                    const double* dbl_buf = mpr.get_dbl_buf();
                    int32_t n_col_est = 10;
                    Eigen::MatrixXd geno_mat(n_overlapping_samples, n_col_est);
                    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> geno_mask(n_overlapping_samples, n_col_est);
                    uint32_t n_geno_missing = 0;
                    int32_t icol = 0;

                    int32_t gcs[3] = {0, 0, 0};
                    do {
                        gcs[0] = gcs[1] = gcs[2] = 0; // initialize genotype counts

                        const plink_var_t& var = mpr.get_current_variant();
                        if ( var.pos > region.end0 ) {
                            break;
                        }
                        std::string cpra_s(var.to_string());
                        mpr.get_genos();

                        if ( icol >= n_col_est ) {
                            geno_mat.conservativeResize(n_overlapping_samples, n_col_est * 2);
                            geno_mask.conservativeResize(n_overlapping_samples, n_col_est * 2);
                            n_col_est *= 2;
                        }
                        int32_t an = 0;
                        double ac = 0;
                        if ( mpr.is_dosage_present() ) {
                            if ( dbl_buf == NULL) {
                                dbl_buf = mpr.get_dbl_buf();
                            }
                            double sumsq = 0;
                            for(int32_t i =0; i < n_overlapping_samples; ++i) {
                                double ds = 2.0 - dbl_buf[i];
                                geno_mat(i, icol) = ds;
                                geno_mask(i, icol) = true; // not missing
                                ac += ds;
                                an += 2;
                                ++gcs[(ds < 0.5) ? 0 : ( (ds < 1.5) ? 1 : 2 )];
                                sumsq += (ds * ds);
                            } 
                            if ( ac == 0 || an == ac ) {
                                // skip monomorphic variants
                                continue;
                            }
                            // E(Var(g)) = 2 * af * (1-af) = 2 * ac * ( an - ac ) / an / an
                            // Var(g) = sumsq / (an / 2) - 4 * ac * ac / an / an
                            // ratio = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) )
                            double info = ( sumsq / 2.0 * an - 4.0 * ac * ac ) / ( 2.0 * ac * ( an - ac ) );
                            infos.push_back(info);
                        }
                        else {
                            for(int32_t i = 0; i < n_overlapping_samples; ++i) {
                                switch(int_buf[i]) { // make sure to convert 1-based index to 0-based
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
                            if ( ac == 0 || an == ac ) {
                                // skip monomorphic variants
                                continue;
                            }
                            double mean = (double)ac / (double)an * 2.0;
                            // Exp(Var(g)) = af * (1-af) * 2 = mean * (2 - mean) / 2;
                            // Var(g) = EX^2 - EX^2 = (4 * n_2 + 1 * n_1)/n - mean^2
                            double info = ((4.0 * gcs[2] + gcs[1])/(an/2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
                            infos.push_back(info);
                            //notice("mean = %.5g, an = %d, ac = %d", mean, an, ac);
                            for(int32_t i = 0; i < n_overlapping_samples; ++i) {
                                switch(int_buf[i]) {
                                case 0:
                                    geno_mat(i, icol) = 2.0 - mean; // homalt
                                    geno_mask(i, icol) = true; // not missing
                                    break;
                                case 1:
                                    geno_mat(i, icol) = 1.0 - mean; // het
                                    geno_mask(i, icol) = true; // not missing
                                    break;
                                case 2:
                                    geno_mat(i, icol) = 0.0 - mean; // homref
                                    geno_mask(i, icol) = true; // not missing
                                    break;
                                default:
                                    geno_mat(i, icol) = 0; // missing - mean imputation
                                    geno_mask(i, icol) = false; // missing
                                    ++n_geno_missing;
                                    break;
                                }
                            }
                        }
                        v_cpra.push_back(cpra_t(cpra_s.c_str())); 
                        acs.push_back(ac);
                        ans.push_back(an);
                        gc0s.push_back(gcs[0]);
                        gc1s.push_back(gcs[1]);
                        gc2s.push_back(gcs[2]);
                        ++icol;

                        if ( icol >= max_chunk_vars ) { // process the current chunk perform association mapping
                            notice("Processing %d variants for phenotype %s", icol, pheno_matrix.pheno_ids[phe_idx].c_str());
                            
                            geno_mat.conservativeResize(n_overlapping_samples, icol);
                            geno_mask.conservativeResize(n_overlapping_samples, icol);

                            assoc_single_trait( 
                                wf, // output file handle
                                pheno_matrix.pheno_ids[phe_idx].c_str(), // phenotype ID
                                phe_vec,      // phenotype vector
                                phe_rint_vec, // rinted phenotype vector 
                                phe_mask_vec, // phenotype mask vector
                                phe_has_missing, // if the phenotype has missing values
                                skip_rint, // skip rank-based inverse normal transformation
                                geno_mat,    // genotype matrix
                                geno_mask,   // genotype mask matrix
                                n_geno_missing > 0, // if the genotype has missing values
                                v_cpra,      // variant pairs
                                ans,         // allele counts
                                acs,         // allele counts
                                gc0s,        // genotype counts for genotype 0
                                gc1s,        // genotype counts for genotype 1
                                gc2s,        // genotype counts for genotype 2
                                infos       // information values
                            );
                            icol = 0;
                            n_col_est = 10;
                            ans.clear();
                            acs.clear();
                            gc0s.clear();
                            gc1s.clear();
                            gc2s.clear();
                            infos.clear();
                            v_cpra.clear();
                            geno_mat.resize(n_overlapping_samples, n_col_est);
                            geno_mask.resize(n_overlapping_samples, n_col_est);
                            n_geno_missing = 0;
                        }
                    }
                    while ( mpr.read_pivar() );

                    if ( icol > 0 ) {
                        if ( icol < n_col_est ) {
                            geno_mat.conservativeResize(n_overlapping_samples, icol);
                            geno_mask.conservativeResize(n_overlapping_samples, icol);
                        }

                        notice("Processing %d variants for phenotype %s", icol, pheno_matrix.pheno_ids[phe_idx].c_str());
                        
                        geno_mat.conservativeResize(n_overlapping_samples, icol);
                        geno_mask.conservativeResize(n_overlapping_samples, icol);

                        assoc_single_trait( 
                            wf, // output file handle
                            pheno_matrix.pheno_ids[phe_idx].c_str(), // phenotype ID
                            phe_vec,      // phenotype vector
                            phe_rint_vec, // rinted phenotype vector 
                            phe_mask_vec, // phenotype mask vector
                            phe_has_missing, // if the phenotype has missing values
                            skip_rint, // skip rank-based inverse normal transformation
                            geno_mat,    // genotype matrix
                            geno_mask,   // genotype mask matrix
                            n_geno_missing > 0, // if the genotype has missing values
                            v_cpra,      // variant pairs
                            ans,         // allele counts
                            acs,         // allele counts
                            gc0s,        // genotype counts for genotype 0 
                            gc1s,        // genotype counts for genotype 1
                            gc2s,        // genotype counts for genotype 2
                            infos       // information values
                        );
                        icol = 0;
                        ans.clear();
                        acs.clear();
                        gc0s.clear();
                        gc1s.clear();
                        gc2s.clear();
                        infos.clear();
                        v_cpra.clear();
                        geno_mat.resize(n_overlapping_samples, n_col_est);
                        geno_mask.resize(n_overlapping_samples, n_col_est);
                        n_geno_missing = 0;
                    }
                }
            }
        }   
    }
    else {
        for(pair_map_it = pair_map.begin(); pair_map_it != pair_map.end(); ++pair_map_it) {
            int32_t phe_idx = pair_map_it->first;
            const std::set<cpra_t>& pairs = pair_map_it->second;
            if ( pairs.empty() ) continue;

            notice("Processing %zu pairs for phenotype %s", pairs.size(), pheno_matrix.pheno_ids[phe_idx].c_str());
            int32_t n_col_est = 10;
            Eigen::MatrixXd geno_mat(n_overlapping_samples, n_col_est);
            Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> geno_mask(n_overlapping_samples, n_col_est);
            uint32_t n_geno_missing = 0;
            int32_t icol = 0;
            Eigen::VectorXd phe_vec(n_overlapping_samples);
            //phe_vec = phe_mat_adj.row(phe_idx); // get the phenotype vector for the current trait
            phe_vec = pheno_matrix.pheno_mat.col(phe_idx);

            //Eigen::Vector<bool, Eigen::Dynamic> phe_mask_vec = phe_mask.row(phe_idx); // get the missing values for the current trait
            Eigen::Vector<bool, Eigen::Dynamic> phe_mask_vec = pheno_matrix.pheno_mask.col(phe_idx); // get the missing values for the current trait

            bool phe_has_missing = !phe_mask_vec.all();
            Eigen::VectorXd phe_rint_vec;
            if ( ! skip_rint ) {
                // perform rank-based inverse normal transformation
                notice("Performing rank-based inverse normal transformation for phenotype %s", pheno_matrix.pheno_ids[phe_idx].c_str());
                if ( phe_has_missing ) {
                    phe_rint_vec = rint_with_missing(phe_vec, pheno_matrix.pheno_mask.col(phe_idx));
                }
                else {
                    phe_rint_vec = rint_without_missing(phe_vec);
                }
                if ( phe_rint_vec.size() != n_overlapping_samples) {
                    error("Rank-based inverse normal transformation failed for phenotype %s", pheno_matrix.pheno_ids[phe_idx].c_str());
                }
            }
            else {
                notice("Skipping rank-based inverse normal transformation for phenotype %s", pheno_matrix.pheno_ids[phe_idx].c_str());
            }

            // for(int32_t i = 0; i < n_row; ++i) {
            //     //phe_vec(i) = phe_mat_adj(phe_idx, phe_idxs[i]);
            //     phe_vec(i) = phe_mat_adj(phe_idx, i);
            // }
            // read the genotypes for each pair
            std::vector<cpra_t> v_cpra;
            std::vector<int32_t> ans;
            std::vector<double> acs;
            std::vector<int32_t> gc0s, gc1s, gc2s;
            std::vector<double> infos;
            const std::vector<int32_t>& int_buf = mpr.get_int_buf();
            const double* dbl_buf = mpr.get_dbl_buf();
            // notice("int_buf.size() = %zu", int_buf.size());
            for(std::set<cpra_t>::iterator it = pairs.begin(); it != pairs.end(); ++it) {
                std::string cpra_s(it->to_string());
                //notice("Reading variant %s", cpra_s.c_str());
                if ( !mpr.read_pivar(cpra_s.c_str()) ) {
                    notice("Skipping pair %s\t%s, which is not found in the pvar file", it->to_string().c_str(), pheno_matrix.pheno_ids[phe_idx].c_str());
                    continue;
                }   
                notice("Reading genotypes for variant %s", cpra_s.c_str());
                mpr.get_genos();
                notice("Finished reading genotypes for variant %s", cpra_s.c_str());
                // construct the input for association analysis
                if ( icol >= n_col_est ) {
                    geno_mat.conservativeResize(n_overlapping_samples, n_col_est * 2);
                    geno_mask.conservativeResize(n_overlapping_samples, n_col_est * 2);
                    n_col_est *= 2;
                }
                int32_t an = 0;
                double ac = 0;
                int32_t gcs[3] = {0, 0, 0};
                if ( mpr.is_dosage_present() ) {
                    if ( dbl_buf == NULL) {
                        dbl_buf = mpr.get_dbl_buf();
                    }
                    double sumsq = 0;
                    for(int32_t i =0; i < n_overlapping_samples; ++i) {
                        //notice("Dosage[%d] = %.5g", i, dbl_buf[i]);
                        double ds = 2.0 - dbl_buf[i];
                        geno_mat(i, icol) = ds;
                        geno_mask(i, icol) = true; // not missing
                        ac += ds;
                        an += 2;
                        ++gcs[(ds < 0.5) ? 0 : ( (ds < 1.5) ? 1 : 2 )];
                        sumsq += (ds * ds);
                    } 
                    if ( ac == 0 || an == ac ) {
                        // skip monomorphic variants
                        continue;
                    }

                    // E(Var(g)) = 2 * af * (1-af) = 2 * ac * ( an - ac ) / an / an
                    // Var(g) = sumsq / (an / 2) - 4 * ac * ac / an / an
                    // ratio = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) )
                    double info = ( sumsq / 2 * an - 4 * ac * ac ) / ( 2 * ac * ( an - ac ) );
                    infos.push_back(info);
                }
                else {
                    //int32_t gcs[3] = {0, 0, 0};
                    for(int32_t i = 0; i < n_overlapping_samples; ++i) {
                        switch(int_buf[i]) { // make sure to convert 1-based index to 0-based
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

                    if ( ac == 0 || an == ac ) {
                        // skip monomorphic variants
                        continue;
                    }

                    double mean = (double)ac / (double)an * 2.0;
                    // Exp(Var(g)) = af * (1-af) * 2 = mean * (2 - mean) / 2;
                    // Var(g) = EX^2 - EX^2 = (4 * n_2 + 1 * n_1)/n - mean^2
                    double info = ((4.0 * gcs[2] + gcs[1])/(an/2.0) - mean * mean) / (mean * (2.0 - mean) / 2.0);
                    infos.push_back(info);
                    //notice("mean = %.5g, an = %d, ac = %d", mean, an, ac);
                    for(int32_t i = 0; i < n_overlapping_samples; ++i) {
                        switch(int_buf[i]) {
                        case 0:
                            geno_mat(i, icol) = 2.0 - mean; // homalt
                            geno_mask(i, icol) = true; // not missing
                            break;
                        case 1:
                            geno_mat(i, icol) = 1.0 - mean; // het
                            geno_mask(i, icol) = true; // not missing
                            break;
                        case 2:
                            geno_mat(i, icol) = 0.0 - mean; // homref
                            geno_mask(i, icol) = true; // not missing
                            break;
                        default:
                            geno_mat(i, icol) = 0; // missing - mean imputation
                            geno_mask(i, icol) = false; // missing
                            ++n_geno_missing;
                            break;
                        }
                    }
                }
                v_cpra.push_back(*it); 
                acs.push_back(ac);
                ans.push_back(an);
                gc0s.push_back(gcs[0]);
                gc1s.push_back(gcs[1]);
                gc2s.push_back(gcs[2]);
                ++icol;
            }

            if ( icol > 0 ) {
                if ( icol < n_col_est ) {
                    geno_mat.conservativeResize(n_overlapping_samples, icol);
                    geno_mask.conservativeResize(n_overlapping_samples, icol);
                }

                assoc_single_trait( 
                    wf, // output file handle
                    pheno_matrix.pheno_ids[phe_idx].c_str(), // phenotype ID
                    phe_vec,      // phenotype vector
                    phe_rint_vec, // rinted phenotype vector 
                    phe_mask_vec, // phenotype mask vector
                    phe_has_missing, // if the phenotype has missing values
                    skip_rint, // skip rank-based inverse normal transformation
                    geno_mat,    // genotype matrix
                    geno_mask,   // genotype mask matrix
                    n_geno_missing > 0, // if the genotype has missing values
                    v_cpra,      // variant pairs
                    ans,         // allele counts
                    acs,         // allele counts
                    gc0s,        // genotype counts for genotype 0
                    gc1s,        // genotype counts for genotype 1
                    gc2s,        // genotype counts for genotype 2
                    infos       // information values
                );
            }

            // std::vector<slr_sumstat_t> sumstats;
            // std::vector<slr_sumstat_t> sumstats_rint;
            // if ( ( n_geno_missing == 0 ) && !phe_has_missing ) {
            //     simple_linear_regression_without_missing(phe_vec, geno_mat, sumstats);
            //     if ( ! skip_rint ) {
            //         simple_linear_regression_without_missing(phe_rint_vec, geno_mat, sumstats_rint);
            //     }
            // }
            // else {
            //     simple_linear_regression_with_missing(phe_vec, phe_mask_vec, geno_mat, geno_mask, sumstats);
            //     if ( ! skip_rint ) {
            //         simple_linear_regression_with_missing(phe_rint_vec, phe_mask_vec, geno_mat, geno_mask, sumstats_rint);
            //     }
            // }

            // // print the results
            // for(int32_t i=0; i < (int32_t)v_cpra.size(); ++i) {
            //     const slr_sumstat_t& ss = sumstats[i];
            //     hprintf(wf, "%s\t%s\t%d\t%s\t%s\t%s\t%.5g\t%.5g\t%d\tADD\t%.6g\t%.6g\t%.6g\t%.6g",
            //         phe_trait_ids[phe_idx].c_str(), // TRAIT
            //         v_cpra[i].chrom.c_str(), // CHROM
            //         v_cpra[i].pos,           // POS
            //         v_cpra[i].to_string().c_str(), // ID
            //         v_cpra[i].ref.c_str(),  // REF
            //         v_cpra[i].alts.c_str(), // ALT
            //         (double)acs[i] / (double)ans[i], // AF
            //         infos[i],  // INFO - placeholder, not calculated
            //         ss.n_obs, // N - number of samples
            //         ss.beta, // BETA
            //         ss.se,   // SE
            //         ss.tstat, // TSTAT
            //         ss.log10p); // LOG10P
            //     if ( ! skip_rint ) {
            //         const slr_sumstat_t& ss_rint = sumstats_rint[i];
            //         hprintf(wf, "\t%.6g\t%.6g\t%.6g\t%.6g", // BETA_RINT, SE_RINT, TSTAT_RINT, LOG10P_RINT
            //             ss_rint.beta, ss_rint.se, ss_rint.tstat, ss_rint.log10p);
            //     }
            //     hprintf(wf, "\n");
            // }
        }
    }
    hts_close(wf); // close the output file
    notice("Analysis finished");
    return 0;
}