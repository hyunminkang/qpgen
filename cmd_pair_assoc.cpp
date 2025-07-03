#include "qgenlib/params.h"
#include "qgenlib/tsv_reader.h"
#include "qgenlib/qgen_utils.h"
#include "qgenlib/phred_helper.h"
#include "qgenlib/hts_utils.h"
#include "assoc_utils.h"
#include "qpgen.h"
#include "Eigen/Dense"
#include <cmath>

int32_t cmd_pair_assoc(int32_t argc, char **argv)
{
    std::string pgenf;
    std::string psamf;
    std::string pivarf;
    std::string phef;
    std::string covf;
    std::string pairf;
    std::string outf;
    int32_t jump_thres_bp = 1000000; 
    
    paramList pl;

    BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Options", NULL)
    LONG_STRING_PARAM("pgen", &pgenf, "Input PLINK2 genotype file")
    LONG_STRING_PARAM("psam", &psamf, "Input PLINK2 sample file")
    LONG_STRING_PARAM("pivar", &pivarf, "Input PLINK2 index pvar file (bgzipped and tabix)")
    LONG_STRING_PARAM("pheno", &phef, "Input phenotype matrix")
    LONG_STRING_PARAM("pairs", &pairf, "Input file containing pairs of variant and phenotype IDs to be tested")
    LONG_STRING_PARAM("cov", &covf, "Input covariate matrix (optional)")

    LONG_PARAM_GROUP("Output options", NULL)
    LONG_STRING_PARAM("out", &outf, "Output prefix")

    LONG_PARAM_GROUP("Auxiliary options", NULL)
    LONG_INT_PARAM("jump-thres-bp", &jump_thres_bp, "Jump threshold in base pairs for the variant index (default: 1000000)")
    END_LONG_PARAMS();

    pl.Add(new longParams("Available Options", longParameters));
    pl.Read(argc, argv);
    pl.Status();

    // open the Plink file
    PgenIdxReader pr;
    if (!pr.prep_pgen(pgenf.c_str(), pivarf.c_str(), psamf.c_str()))
        error("Failed to prepare pgen files with the following files:\n%s\n%s\n%s", pgenf.c_str(), pivarf.c_str(), psamf.c_str());

    // load the sample IDs first
    if (!pr.load_psam(psamf.c_str()))
        error("Failed to load the sample IDs from %s", psamf.c_str());


    // read the phenotype file
    tsv_reader tr_phe(phef.c_str());
    // read the header line
    if (!tr_phe.read_line())
        error("Cannot read the header line from the phenotype file %s", phef.c_str());
    std::vector<std::string> phe_sample_ids;
    std::map<std::string, int32_t> phe_id2idx;
    for(int32_t i = 4; i < tr_phe.nfields; ++i) {
        phe_sample_ids.push_back(tr_phe.str_field_at(i));
        phe_id2idx[tr_phe.str_field_at(i)] = i - 4; // store the index of the sample ID
    }

    // read the covariate file if exists
    std::set<std::string> overlapping_sample_ids;
    std::vector<int32_t> cov_idxs;
    std::vector<int32_t> phe_idxs;
    std::vector<std::string> cov_ids;

    Eigen::MatrixXd cov_mat;
    if (!covf.empty())
    {
        tsv_reader tr_cov(covf.c_str());
        // read the header line
        if (!tr_cov.read_line())
            error("Cannot read the header line from the covariate file %s", covf.c_str());
        std::map<std::string, int32_t> cov_id2idx;
        for(int32_t i = 1; i < tr_cov.nfields; ++i) {
            const char* id = tr_cov.str_field_at(i);
            // check if the sample ID exists in both
            if ( pr.samp2idx.find(id) != pr.samp2idx.end() ) {   // exists in psam
                if ( phe_id2idx.find(id) != phe_id2idx.end() ) { // exists in pheno
                    overlapping_sample_ids.insert(id);                    
                    cov_id2idx[id] = i - 1; // store the index of the sample ID
                }
            }
        }

        // construct the indices
        std::set<std::string>::iterator it;
        pr.samp_idx.clear();
        for(it = overlapping_sample_ids.begin(); it != overlapping_sample_ids.end(); ++it) {
            cov_idxs.push_back(cov_id2idx[*it]);
            phe_idxs.push_back(phe_id2idx[*it]);
            pr.samp_idx.push_back(pr.samp2idx[*it]+1); 
        }

        notice("Identified %zu overlapping sample IDs between genotype, phenotype and covariate files", overlapping_sample_ids.size());

        // load the covariate data
        int32_t nrow_est = 100;
        cov_mat.resize(nrow_est, cov_idxs.size());
        int32_t irow = 0;
        while( tr_cov.read_line() ) {
            cov_ids.push_back(tr_cov.str_field_at(0)); // first field is the covariate name
            if ( irow >= nrow_est ) {
                cov_mat.conservativeResize(nrow_est * 2, cov_idxs.size());
                nrow_est *= 2;
            }
            for(int32_t i = 0; i < cov_idxs.size(); ++i) {
                cov_mat(irow, i) = tr_cov.double_field_at(cov_idxs[i]+1);
            }
            ++irow;
        }
        cov_mat.conservativeResize(irow, cov_idxs.size());

        notice("Loaded %d covariates for %zu samples", irow, (int32_t)cov_idxs.size());
    }
    else {
        std::map<std::string, int32_t>::iterator it;
        pr.samp_idx.clear();
        for(it = phe_id2idx.begin(); it != phe_id2idx.end(); ++it) {
            if ( pr.samp2idx.find(it->first) != pr.samp2idx.end() ) { // exists in psam
                overlapping_sample_ids.insert(it->first);
                phe_idxs.push_back(it->second);
                pr.samp_idx.push_back(pr.samp2idx[it->first]+1);
            }
        }
        notice("Identified %zu overlapping sample IDs between genotype and phenotype files", overlapping_sample_ids.size());
        notice("phe_idxs.size() = %zu, pr.samp_idx.size() = %zu", phe_idxs.size(), pr.samp_idx.size());
    }

    // for(int32_t i=0; i < 10; ++i) {
    //     notice("%d\t%d\t%d\t%d\t%s\t%s", 
    //         i, phe_idxs[i], cov_idxs.size() > i ? cov_idxs[i] : -1, pr.samp_idx[i],
    //         phe_sample_ids[phe_idxs[i]].c_str(),
    //         pr.samps[pr.samp_idx[i]-1].indID.c_str());
    // }

    // read the phenotype matrix
    std::vector<std::string> phe_chroms;
    std::vector<int32_t> phe_begs;
    std::vector<int32_t> phe_ends;
    std::vector<std::string> phe_ids;
    std::map<std::string, int32_t> phe2idx; // map phenotype ID to idx
    std::map<std::string, int32_t>::iterator phe2idx_it;
    int32_t nrow_est = 100;
    int32_t irow = 0;
    Eigen::MatrixXd phe_mat;
    phe_mat.resize(nrow_est, phe_idxs.size());
    while( tr_phe.read_line() ) {
        phe_chroms.push_back(tr_phe.str_field_at(0));
        phe_begs.push_back(tr_phe.int_field_at(1));
        phe_ends.push_back(tr_phe.int_field_at(2));
        phe_ids.push_back(tr_phe.str_field_at(3)); // the phenotype ID
        phe2idx[tr_phe.str_field_at(3)] = irow; // store the index of the phenotype ID
        if ( irow >= nrow_est ) {
            phe_mat.conservativeResize(nrow_est * 2, phe_idxs.size());
            nrow_est *= 2;
        }
        for(int32_t i = 0; i < phe_idxs.size(); ++i) {
            phe_mat(irow, i) = tr_phe.double_field_at(phe_idxs[i]+4);
        }
        ++irow;
    }
    phe_mat.conservativeResize(irow, phe_idxs.size());
    notice("Loaded %d phenotypes for %zu samples. Dimension : %d x %d", irow, (int32_t)phe_idxs.size(), phe_mat.rows(), phe_mat.cols());

    Eigen::MatrixXd phe_mat_adj;
    if ( covf.empty() ) {
        notice("No covariates provided. Centering the phenotypes by sample means");
        // center the phenotypes by sample means
        phe_mat_adj = phe_mat;
        center_rows(phe_mat_adj);
        //notice("%.5g", phe_mat_adj(0,0));
    }
    else {
        notice("Adjusting the phenotypes by covariates");
        phe_mat_adj = bulk_adjust_for_covariates(phe_mat, cov_mat);
    }
    notice("Finished adjusting the phenotype matrix");

    // read the variant-trait pair lists
    tsv_reader tr_pair(pairf.c_str());
    std::map<int32_t, std::set<cpra_t> > pair_map;
    std::map<int32_t, std::set<cpra_t> >::iterator pair_map_it;
    notice("Reading the variant-trait pairs from %s", pairf.c_str());
    int32_t n_pairs = 0;

    // open the output file gz or plain based on the extension
    htsFile* wf = hts_open(outf.c_str(), outf.substr(outf.length() - 3).compare(".gz") == 0 ? "wz" : "w");
    if ( wf == NULL ) {
        error("Cannot open output file %s for writing", outf.c_str());
    }
    // write the header line
    hprintf(wf, "#TRAIT\tCHROM\tPOS\tID\tREF\tALT\tAF\tINFO\tN\tTEST\tBETA\tSE\tTSTAT\tLOG10P\tEXTRA\n");

    while( tr_pair.read_line() ) {
        const char* var_id = tr_pair.str_field_at(0);
        const char* phe_id = tr_pair.str_field_at(1);
        phe2idx_it = phe2idx.find(phe_id);
        if ( phe2idx_it == phe2idx.end() ) {
            notice("Skipping phenotype %s, which is not observed in %s", phe_id, phef.c_str());
        }
        else {
            int32_t phe_idx = phe2idx_it->second;
            cpra_t cpra(var_id);
            if ( pair_map[phe_idx].insert(cpra).second ) {
                ++n_pairs;
            }
            else {
                notice("Skipping duplicate pair %s\t%s", var_id, phe_id);
            }
        }
    }
    notice("Read %d variant-trait pairs across %zu traits", n_pairs, pair_map.size());

    for(pair_map_it = pair_map.begin(); pair_map_it != pair_map.end(); ++pair_map_it) {
        int32_t phe_idx = pair_map_it->first;
        const std::set<cpra_t>& pairs = pair_map_it->second;
        if ( pairs.empty() ) continue;

        notice("Processing %zu pairs for phenotype %s", pairs.size(), phe_ids[phe_idx].c_str());
        int32_t n_col_est = 10;
        int32_t n_row = (int32_t)overlapping_sample_ids.size();
        Eigen::MatrixXd geno_mat(n_row, n_col_est);
        int32_t icol = 0;
        Eigen::VectorXd phe_vec(n_row);
        for(int32_t i = 0; i < n_row; ++i) {
            phe_vec(i) = phe_mat_adj(phe_idx, phe_idxs[i]);
        }
        // read the genotypes for each pair
        std::vector<cpra_t> v_cpra;
        std::vector<int32_t> ans;
        std::vector<int32_t> acs;
        for(std::set<cpra_t>::iterator it = pairs.begin(); it != pairs.end(); ++it) {
            std::string cpra_s(it->to_string());
            //notice("Reading variant %s", cpra_s.c_str());
            if ( !pr.read_pivar(cpra_s.c_str()) ) {
                notice("Skipping pair %s\t%s, which is not found in the pvar file", it->to_string().c_str(), phe_ids[phe_idx].c_str());
                continue;
            }   
            pr.get_genos();
            // construct the input for association analysis
            if ( icol >= n_col_est ) {
                geno_mat.conservativeResize(n_row, n_col_est * 2);
                n_col_est *= 2;
            }
            int32_t an = 0, ac = 0;
            for(int32_t i = 0; i < n_row; ++i) {
                switch(pr.int_buf[pr.samp_idx[i]-1]) { // make sure to convert 1-based index to 0-based
                case 0:
                    an += 2;
                    ac += 2;
                    break;
                case 1:
                    an += 2;
                    ++ac;
                    break;
                case 2:
                    an += 2;
                    break;
                }
            }
            double mean = (double)ac / (double)an * 2.0;
            //notice("mean = %.5g, an = %d, ac = %d", mean, an, ac);
            for(int32_t i = 0; i < n_row; ++i) {
                switch(pr.int_buf[pr.samp_idx[i]-1]) {
                case 0:
                    geno_mat(i, icol) = 2.0 - mean; // homalt
                    break;
                case 1:
                    geno_mat(i, icol) = 1.0 - mean; // het
                    break;
                case 2:
                    geno_mat(i, icol) = 0.0 - mean; // homref
                    break;
                default:
                    geno_mat(i, icol) = 0; // missing - mean imputation
                    break;
                }
            }
            v_cpra.push_back(*it); 
            acs.push_back(ac);
            ans.push_back(an);
            ++icol;
        }
        if ( icol < n_col_est ) {
            geno_mat.conservativeResize(n_row, icol);
        }

        // for(int32_t i=0; i < n_row; ++i) {
        //     printf("%s\t%.5g\t%.5g\n", phe_sample_ids[phe_idxs[i]].c_str(), phe_vec(i), geno_mat(i, 0));
        // }
        
        //std::cout << geno_mat << std::endl;
        // Perform linear regression for each pair
        std::vector<slr_sumstat_t> sumstats;
        simple_linear_regression(phe_vec, geno_mat, sumstats);

        // print the results
        for(int32_t i=0; i < (int32_t)v_cpra.size(); ++i) {
            const slr_sumstat_t& ss = sumstats[i];
            hprintf(wf, "%s\t%s\t%d\t%s\t%s\t%s\t%.5g\t%.5g\t%d\tADD\t%.6g\t%.6g\t%.6g\t%.6g\tNA\n",
                phe_ids[phe_idx].c_str(), // TRAIT
                v_cpra[i].chrom.c_str(), // CHROM
                v_cpra[i].pos,           // POS
                v_cpra[i].to_string().c_str(), // ID
                v_cpra[i].ref.c_str(),  // REF
                v_cpra[i].alts.c_str(), // ALT
                (double)acs[i] / (double)ans[i], // AF
                1.000,  // INFO - placeholder, not calculated
                ans[i]/2, // N - number of samples
                ss.beta, // BETA
                ss.se,   // SE
                ss.tstat, // TSTAT
                ss.log10p); // LOG10P
            // notice("%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\n",
            //         cpra_strs[i].c_str(), phe_ids[phe_idx].c_str(),
            //         ss.beta, ss.se, ss.tstat, ss.log10p);
        }
    }
    hts_close(wf); // close the output file
    notice("Analysis finished");
    return 0;
}
