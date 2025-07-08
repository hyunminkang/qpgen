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
    //int32_t ichrom;
    std::string schrom;
    int32_t pos;
    std::string vid;
    std::string ref;
    std::vector<std::string> alts;
    double cM;     // centimorgan

    std::string to_string(char delim = ':') const {
        std::string s;
        s.append(schrom);
        s.push_back(delim);
        s.append(std::to_string(pos));
        s.push_back(delim);
        s.append(ref);
        s.push_back(delim);
        for(int32_t i = 0; i < (int32_t)alts.size(); ++i) {
            if ( i > 0 ) s.push_back(',');
            s.append(alts[i]);
        }
        return s;
    }

    // additional information that can be filled with genotype data
    uint32_t an;    // number of non-missing alleles
    uint32_t ns;    // number of non-missing genotypes
    std::vector<uint32_t> acs; // allele counts (including ref)
    std::vector<uint32_t> gcs; // genotype counts
    std::vector<double> afs;
};
typedef struct _plink_var_t plink_var_t;

struct _cbe_t {
    std::string chrom;
    int32_t beg1;
    int32_t end0; // 1-based inclusive start and 0-based exclusive end

    bool operator<(const struct _cbe_t& other) const {
        if ( chrom < other.chrom ) return true;
        if ( chrom > other.chrom ) return false;
        if ( beg1 < other.beg1 ) return true;
        if ( beg1 > other.beg1 ) return false;
        return end0 < other.end0; 
    }

    std::string to_string(char delim1 = ':', char delim2 = '-') const {
        std::string s;
        s.append(chrom);
        s.push_back(delim1);
        s.append(std::to_string(beg1));
        s.push_back(delim2);
        s.append(std::to_string(end0));
        return s;
    }

    _cbe_t(const char* _chrom, int32_t _beg1, int32_t _end0) 
        : chrom(_chrom), beg1(_beg1), end0(_end0) {
        if ( _beg1 < 1 || _end0 < 0 || _end0 <= _beg1 ) {
            error("Invalid CBE format: %s:%d-%d", _chrom, _beg1, _end0);
        }
    }

    _cbe_t(const char* s) {
        const char delim1 = ':';
        const char delim2 = '-';
        const char* p = strchr(s, delim1);
        if ( p == NULL ) {
            error("Invalid CBE format: %s", s);
        }
        chrom.assign(s, p - s);
        s = p + 1;
        p = strchr(s, delim2);
        if ( p == NULL ) {
            error("Invalid CBE format: %s", s);
        }
        beg1 = atoi(s);
        s = p + 1;
        if ( *s == '\0' ) {
            error("Invalid CBE format: %s", s);
        }
        else {
            end0 = atoi(s); // 0-based exclusive end
        }
    }
};

typedef struct _cbe_t cbe_t;

struct _cpra_t {
    std::string chrom;
    int32_t pos;
    std::string ref;
    std::string alts;

    bool operator<(const struct _cpra_t& other) const {
        if ( chrom < other.chrom ) return true;
        if ( chrom > other.chrom ) return false;
        if ( pos < other.pos ) return true;
        if ( pos > other.pos ) return false;
        if ( ref < other.ref ) return true;
        if ( ref > other.ref ) return false;
        return alts < other.alts; // compare alts last
    }

    std::string to_string(char delim = ':') const {
        std::string s;
        s.append(chrom);
        s.push_back(delim);
        s.append(std::to_string(pos));
        s.push_back(delim);
        s.append(ref);
        s.push_back(delim);
        s.append(alts);
        return s;
    }

    _cpra_t(const char* s, char delim = ':') {
        const char* p = strchr(s, delim);
        if ( p == NULL ) {
            error("Invalid CPRA format: %s", s);
        }
        chrom.assign(s, p - s);
        s = p + 1;
        p = strchr(s, delim);
        if ( p == NULL ) {
            error("Invalid CPRA format: %s", s);
        }
        pos = atoi(s);
        s = p + 1;
        p = strchr(s, delim);
        if ( p == NULL ) {
            error("Invalid CPRA format: %s", s);
        }
        ref.assign(s, p - s);
        s = p + 1;
        if ( *s == '\0' ) {
            error("Invalid CPRA format: %s", s);
        }
        else {
            alts.assign(s); // the rest is alt alleles
        }
        // while ( ( p = strchr(s, ',') ) != NULL ) {
        //     alts.push_back(std::string(s, p - s));
        //     s = p + 1;
        // }
        // if ( *s != '\0' ) { // last alt
        //     alts.push_back(s);
        // }
    }
};
typedef struct _cpra_t cpra_t;

struct _plink_samp_t {
    std::string famID;
    std::string indID;
    std::string dadID;
    std::string momID;
    int32_t sex;
    double pheno; // if phenotype is available
};

typedef struct _plink_samp_t plink_samp_t;

// PLINK 2.0 reader with indexed pvar
class PgenIdxReader {
protected:
    PgenReader pgr;
    int32_t nthreads;
    tsv_reader tr_pivar; // read tabixed/indexed pvar file
    double* dbl_buf;
    std::vector<int32_t> int_buf;

    std::string pgenf;
    std::string psamf;
    std::string pivarf;

    // internal variables, should not be modified after calling get_genos()
    std::vector<plink_samp_t> samps;           // sample IDs -- fully loaded in memory 
    std::map<std::string, int32_t> samp2idx;  // sample ID maps 1-based index
    std::vector<int32_t> samp_idx;             // sample indices to subset and load - 1-based index

    int jump_thres_bp;

    bool pivar_loaded;
    bool pgen_loaded;
    bool dosage_present; // true if dosage is present in the pgen file
    plink_var_t cur_var;  // current variant information
    int32_t cur_var_idx;  // current variant index
    int32_t icol_pivar_idx; // column index for the variant ID in the pvar file (0-based)

public:
    PgenIdxReader() : pivar_loaded(false), pgen_loaded(false), nthreads(1), dbl_buf(NULL), cur_var_idx(-1), icol_pivar_idx(8), jump_thres_bp(10000), dosage_present(false) {}
    ~PgenIdxReader() {
        if ( dbl_buf ) free(dbl_buf);
    }

    // functions to load ALL sample and FIRST variant info
    bool prep_pgen(const char* _pgenf, const char* _pivarf, const char* _psamf);

    // functions to set filters for samples and variants - should run after prep
    //void set_filter_sample_id(std::vector<std::string>& samp_ids, bool exclude = false);
    void subset_sample_ids(const std::vector<std::string>& samp_ids, bool exclude = false);
    void subset_sample_indices(const std::vector<int32_t>& samp_indices, bool exclude = false);

    bool read_pos(const char* chrom, int32_t pos); // change the current variant position to a specific CPRA
    bool read_pivar(const char* cpra = NULL); // change the current variant position to a specific CPRA
    bool get_genos(int32_t var_idx = -1);   // read the genotypes at the current variant position
    bool load_psam(const char* _psamf);
    
    int32_t get_n_threads() const { return nthreads; }
    void set_n_threads(int32_t n) { nthreads = n; }

    int32_t get_jump_thres_bp() const { return jump_thres_bp; }
    void set_jump_thres_bp(int32_t thres) { jump_thres_bp = thres; }

    int32_t get_icol_pivar_idx() const { return icol_pivar_idx; }
    void set_icol_pivar_idx(int32_t idx) { icol_pivar_idx = idx; }

    const std::vector<int32_t>& get_int_buf() const { return int_buf; }
    const double* get_dbl_buf() const { return dbl_buf; }
    const std::map<std::string, int32_t>& get_samp2idx() { return samp2idx; }

    bool is_pivar_loaded() const { return pivar_loaded; }
    bool is_pgen_loaded() const { return pgen_loaded; }
    bool is_dosage_present() const { return dosage_present; }
    const plink_var_t& get_current_variant() const { return cur_var; }
    int32_t get_current_variant_idx() const { return cur_var_idx; }

    int32_t get_sample_count() const { return (int32_t)samps.size(); }
    const std::vector<plink_samp_t>& get_samples() const { return samps; }
};

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
    std::map<std::string, int32_t> samp2idx;    // sample ID maps
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