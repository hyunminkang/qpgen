#ifndef __ASSOC_PHENO_H
#define __ASSOC_PHENO_H

#include <set>
#include <string>

#include "qgenlib/qgen_error.h"
#include "qgenlib/hts_utils.h"
#include "qgenlib/genome_loci.h"
#include "qpgen.h"
#include "Eigen/Dense"

class PhenoMatrix {
public:
    Eigen::MatrixXd pheno_mat;
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> pheno_mask;
    std::vector<std::string> pheno_ids;
    std::vector<genomeLocus> pheno_loci;
    genomeLocusMap<int32_t> pheno_locusmap;
    std::vector<std::string> samp_ids;
    std::string format;
    bool has_missing;
    bool has_loci;
    std::set<std::string> missing_strs; 
    std::map<std::string, int32_t> pheno_id2idx;
    std::map<std::string, int32_t> samp_id2idx;

    PhenoMatrix() : has_missing(false), has_loci(false) {}
    ~PhenoMatrix() {}
    bool is_missing(const std::string& str) { return missing_strs.find(str) != missing_strs.end(); }
    bool is_missing(const char* str) { return is_missing(std::string(str)); }
    bool add_missing_str(const std::string& str) { return missing_strs.insert(str).second; }
    void reset_missing_strs() { missing_strs.clear(); }
    
    bool load_pheno_matrix(const char* pheno_file, const char* pheno_format, const char delim = '\0');
    int32_t subset_sample_ids(const std::vector<std::string>& samp_ids);
    int32_t subset_pheno_ids(const std::vector<std::string>& pheno_ids);
    bool subset_sample_pheno_ids(const std::vector<std::string>& samp_ids, const std::vector<std::string>& pheno_ids);
    int32_t subset_sample_indices(const std::vector<int32_t>& samp_indices);
    int32_t subset_pheno_indices(const std::vector<int32_t>& pheno_indices);
    bool subset_sample_pheno_indices(const std::vector<int32_t>& samp_indices, const std::vector<int32_t>& pheno_indices);

    bool sample_ids_sorted() const;
    bool pheno_ids_sorted() const;
};

#endif // __ASSOC_PHENO_H