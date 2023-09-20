#ifndef __QPGEN_H
#define __QPGEN_H

#include "qgenlib/qgen_error.h"
#include "qgenlib/tsv_reader.h"
#include "pgenlibr.h"

// A Reader for PLINK 1.9/2.0 .bed/.pgen file
// with support for tabixed variant index
// NOTE: Many functions are modified versions from regenie at
// https://github.com/rgcgithub/regenie/tree/master/external_libs/pgenlibr
// shared under GNU Lesser General Public License v3.0

struct _plink_var_t {
    int32_t ichrom;
    std::string schrom;
    int32_t pos;
    std::string vid;
    std::string ref;
    std::vector<std::string> alts;
    double cM;     // centimorgan

    // additional information that can be filled with genotype data
    uint32_t an;    // number of non-missing alleles
    uint32_t ns;    // number of non-missing genotypes
    std::vector<uint32_t> acs; // allele counts (including ref)
    std::vector<uint32_t> gcs; // genotype counts
    std::vector<double> afs;
};
typedef struct _plink_var_t plink_var_t;

struct _plink_samp_t {
    std::string famID;
    std::string indID;
    std::string dadID;
    std::string momID;
    int32_t sex;
    double pheno; // if phenotype is available
};

typedef struct _plink_samp_t plink_samp_t;

class PlinkReader {
public:
    bool mode_bed; // true : bed mode, false : pgen mode
    PgenReader pgr; // stores genotypes
    bool pgr_loaded; // true : pgr is loaded, false : pgr is empty
    int32_t nthreads;
    tsv_reader tr_var; // stores variant info
    double* dbl_buf;
    std::vector<int32_t> int_buf;

    // variables for pgen mode
    std::string pgenf;
    std::string psamf;
    std::string pvarf;

    // variables for bed mode
    std::string bedf;
    std::string bimf;
    std::string famf;

    // internal variables
    std::vector<plink_samp_t> samps;             // sample IDs -- fully loaded in memory 
    std::map<std::string, uint32_t> samp2idx;    // sample ID maps
    std::vector<int32_t> samp_idx;              // sample indices to load 

    uint32_t nvars;                     // variant index offset (if removed from memory)
    std::vector<plink_var_t> vars;               // variant IDs -- partially loaded in memory
    //std::map<std::string, std::vector<uint32_t> > var2idx;   // variant ID maps
    //std::map<std::string, uint32_t> cpra2idx;    // CPRA ID maps    

    PlinkReader() : mode_bed(false), pgr_loaded(false), nvars(0), nthreads(1), dbl_buf(NULL) {}

    // functions to load ALL sample and FIRST variant info
    bool prep_prefix(const char* _prefix);
    bool prep_pgen(const char* _pgenf, const char* _pvarf, const char* _psamf);
    bool prep_bed(const char* _bedf, const char* _bimf, const char* _famf);

    // functions to set filters for samples and variants - should run after prep
    void set_filter_sample_id(std::vector<std::string>& samp_ids, bool exclude = false);
    //bool filter_var_cpra(std::vector<std::string>& var_cpra, bool exclude = false);

    // move the variant position to somewhere else
    //bool move_var_idx(uint32_t new_var_idx);
    //bool move_var_pos(const char* chrom, uint32_t pos);

    int32_t get_variant_idx_from_cpra(const std::vector<std::string>& var_cpra, std::vector<int32_t>& variant_idx);
    int32_t pvar_cpra2idx(std::map<std::string, int32_t>& cpra2idx, std::vector<int32_t>& variant_idx);
    int32_t bim_cpra2idx(std::map<std::string, int32_t>& cpra2idx, std::vector<int32_t>& variant_idx);

    bool get_bed_genos_at(int32_t var_idx);
    bool get_pgen_genos_at(int32_t var_idx);
    bool get_genos_at(int32_t var_idx) { if (mode_bed) return get_bed_genos_at(var_idx); else return get_pgen_genos_at(var_idx); }

    // function to load genotypes - should be done after prep and filter
    bool stream_pgen_genos(); // read one marker at a time
    bool stream_bed_genos();
    bool stream_genos() { if (mode_bed) return stream_bed_genos(); else return stream_pgen_genos(); }

    bool load_psam(const char* _psamf);
    bool load_fam(const char* _famf);
    bool stream_pvar();
    bool stream_bim();
};

#endif // __QPGEN_H