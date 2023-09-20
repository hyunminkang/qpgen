#include "qpgen.h"
#include "qgenlib/qgen_utils.h"
#include <cstring>

// load PLINK files in a flexible manner.
// First, it searches for .pgen/.psam/.pvar files
// If not found, it searches for .bed/.bim/.fam files
bool PlinkReader::prep_prefix(const char*  _prefix) {
    // check if .pgen/.psam/.pvar files exist
    std::string str_prefix(_prefix);
    std::string _pgenf = str_prefix + ".pgen";
    std::string _psamf = str_prefix + ".psam";
    std::string _pvarf = str_prefix + ".pvar";

    if (check_file_existence(_pgenf.c_str()) && check_file_existence(_pvarf.c_str()) && check_file_existence(_psamf.c_str())) {
        return prep_pgen(_pgenf.c_str(), _pvarf.c_str(), _psamf.c_str());
    }
    // check if .bed/.bim/.fam files exist
    std::string _bedf = str_prefix + ".bed";
    std::string _bimf = str_prefix + ".bim";
    std::string _famf = str_prefix + ".fam";
    if (check_file_existence(_bedf.c_str()) && check_file_existence(_bimf.c_str()) && check_file_existence(_famf.c_str())) {
        return prep_bed(_bedf.c_str(), _bimf.c_str(), _famf.c_str());
    }
    // if none of the above, return false
    return false;
}

bool PlinkReader::prep_pgen(const char* _pgenf, const char* _pvarf, const char* _psamf) {
    // read the sample file first
    if ( !load_psam(_psamf) ) {
        error("Cannot parse the sample info at %s", _psamf);
    }

    mode_bed = false;   
    pgenf.assign(_pgenf);
    psamf.assign(_psamf);
    pvarf.assign(_pvarf);

    // load the variant file, but only partially
    //if ( !read_pvar(_pvarf) ) {
    //    error("Cannot parse the first variant at %s", _pvarf);
    //}

    // DO NOT LOAD the genotype file yet    
    return true;
}

bool PlinkReader::prep_bed(const char* _bedf, const char* _bimf, const char* _famf) {
    // read the sample file first
    if ( !load_fam(_famf) ) {
        error("Cannot parse the sample info at %s", _famf);
    }

    mode_bed = true;
    bedf.assign(_bedf);
    bimf.assign(_bimf);
    famf.assign(_famf);

    // load the variant file, but only partially
    //if ( !read_bim(_bimf) ) {
    //    error("Cannot parse the first variant at %s", _bimf);
    //}

    // DO NOT LOAD the genotype file yet    
    return true;
}

bool PlinkReader::load_psam(const char* _psamf) {
    tsv_reader tr_ind;

    if ( !tr_ind.open(_psamf) ) { // return false if file cannot be opened
        return false;
    }

    bool has_header = false;
    uint32_t nsamps = 0;
    int32_t idx_fid = -1, idx_iid = -1, idx_pat = -1, idx_mat = -1, idx_sex = -1, idx_pheno = -1;
    while( tr_ind.read_line() ) {
        if ( ( nsamps == 0 ) && ( !has_header ) ) { // first line to process
            if ( tr_ind.str_field_at(0)[0] == '#' ) {  // If header column exists
                for(int32_t i=0; i < tr_ind.nfields; ++i) {
                    const char* str = tr_ind.str_field_at(i) + ( i == 0 ? 1 : 0 );
                    if ( strcmp("FID", str) == 0 ) { idx_fid = i; }
                    else if ( strcmp("IID", str) == 0 ) { idx_iid = i; }
                    else if ( strcmp("PAT", str) == 0 ) { idx_pat = i; }
                    else if ( strcmp("MAT", str) == 0 ) { idx_mat = i; }
                    else if ( strcmp("SEX", str) == 0 ) { idx_sex = i; }
                    else if ( strcmp("PHENO", str) == 0 ) { idx_pheno = i; }
                }
                has_header = true;

                //notice("%d %d %d %d %d %d", idx_fid, idx_iid, idx_pat, idx_mat, idx_sex, idx_pheno);
                continue;
            }
            else {  // if header column does not exist, use the default order
                idx_fid = tr_ind.nfields > 0 ? 0 : -1;
                idx_iid = tr_ind.nfields > 1 ? 1 : -1;
                idx_pat = tr_ind.nfields > 2 ? 2 : -1;
                idx_mat = tr_ind.nfields > 3 ? 3 : -1;
                idx_sex = tr_ind.nfields > 4 ? 4 : -1;
                idx_pheno = tr_ind.nfields > 5 ? 5 : -1;
            }
        }

        // create a sample 
        plink_samp_t samp;
        if ( idx_fid >= 0 ) { samp.famID.assign(tr_ind.str_field_at(idx_fid)); } 
        if ( idx_iid >= 0 ) { samp.indID.assign(tr_ind.str_field_at(idx_iid)); }
        if ( idx_pat >= 0 ) { samp.dadID.assign(tr_ind.str_field_at(idx_pat)); }
        if ( idx_mat >= 0 ) { samp.momID.assign(tr_ind.str_field_at(idx_mat)); }
        if ( idx_sex >= 0 ) { samp.sex = tr_ind.int_field_at(idx_sex); }
        if ( idx_pheno >= 0 ) { samp.pheno = tr_ind.double_field_at(idx_pheno); } 

        // add the sample to the list
        samps.push_back(samp);
        std::string iid = ( samp.famID.empty() || samp.famID.compare("0") == 0 || samp.famID == samp.indID ) ? samp.indID : (samp.famID + "_" + samp.indID);
        samp2idx[iid] = nsamps;
        //error("%s %s %s %d %d", iid.c_str(), samp.famID.c_str(), samp.indID.c_str(), idx_fid, idx_sex);
        ++nsamps;
    }

    // fill in the sample indices to load
    for(uint32_t i=0; i < nsamps; ++i) {
        samp_idx.push_back((int32_t)(i+1));
    }

    return nsamps > 0;
}

bool PlinkReader::load_fam(const char* _famf) {
    tsv_reader tr_ind;

    if ( !tr_ind.open(_famf) ) { // return false if file cannot be opened
        return false;
    }

    uint32_t nsamps = 0;
    while( tr_ind.read_line() ) {
        // create a sample 
        plink_samp_t samp;

        samp.famID.assign(tr_ind.str_field_at(0));
        samp.indID.assign(tr_ind.str_field_at(1));
        if ( tr_ind.nfields > 2 ) { samp.dadID.assign(tr_ind.str_field_at(2)); }
        if ( tr_ind.nfields > 3 ) { samp.momID.assign(tr_ind.str_field_at(3)); }
        if ( tr_ind.nfields > 4 ) { samp.sex = tr_ind.int_field_at(4); }
        if ( tr_ind.nfields > 5 ) { samp.pheno = tr_ind.double_field_at(5); }

        // add the sample to the list
        samps.push_back(samp);
        std::string iid = ( samp.famID.empty() || samp.famID.compare("0") == 0 || samp.famID == samp.indID ) ? samp.indID : (samp.famID + "_" + samp.indID);
        //std::string iid = samp.famID.empty() ? samp.indID : samp.famID + "_" + samp.indID;
        samp2idx[iid] = nsamps;
        ++nsamps;
    }

    // fill in the sample indices to load
    for(uint32_t i=0; i < nsamps; ++i) {
        samp_idx.push_back((int32_t)(i+1));
    }
    return true;
}

void PlinkReader::set_filter_sample_id(std::vector<std::string>& samp_ids, bool exclude) {
    samp_idx.clear();

    std::set<uint32_t> idxset;
    for(int32_t i=0; i < samp_ids.size(); ++i) {
        std::map<std::string, uint32_t>::iterator it = samp2idx.find(samp_ids[i]);
        if ( it != samp2idx.end() ) {
            idxset.insert(it->second);
        }
    }

    for(uint32_t i=0; i < samps.size(); ++i) {
        if ( exclude ) {
            if ( idxset.find(i) == idxset.end() ) {
                samp_idx.push_back(i+1);
            }
        }
        else {
            if ( idxset.find(i) != idxset.end() ) {
                samp_idx.push_back(i+1);
            }
        }
    }

    if ( samp_idx.empty() ) {
        error("No samples to be included after subsetting to %zu", samp_ids.size());
    }
}

bool PlinkReader::stream_pvar() {
    if ( nvars == 0 ) {
        if ( !tr_var.open(pvarf.c_str()) ) { // return false if file cannot be opened
            error("Cannot open %s", pvarf.c_str());
            return false;
        } 
    }

    // ignore headers
    int32_t idx_chrom = -1, idx_pos = -1, idx_id = -1, idx_ref = -1, idx_alt = -1, idx_qual = -1, idx_cm = -1, idx_filter = -1, idx_format = -1, idx_info = -1;
    while( tr_var.read_line() ) {
        const char* s = tr_var.str_field_at(0);
        if ( s[0] != '#' ) {
            break;
        }
        else if ( strncmp(s, "#CHROM", 6) != 0 ) {  // meta line
            continue;
        }

        // this is a header line that contains column names
        for(int32_t i=0; i < tr_var.nfields; ++i) {
            const char* str = tr_var.str_field_at(i) + ( i == 0 ? 1 : 0 );
            if ( strcmp("CHROM", str) == 0 ) { idx_chrom = i; }
            else if ( strcmp("POS", str) == 0 ) { idx_pos = i; }
            else if ( strcmp("ID", str) == 0 ) { idx_id = i; }
            else if ( strcmp("REF", str) == 0 ) { idx_ref = i; }
            else if ( strcmp("ALT", str) == 0 ) { idx_alt = i; }
            else if ( strcmp("QUAL", str) == 0 ) { idx_qual = i; }
            else if ( strcmp("FILTER", str) == 0 ) { idx_filter = i; }
            else if ( strcmp("INFO", str) == 0 ) { idx_info = i; }
            else if ( strcmp("FORMAT", str) == 0 ) { idx_format = i; }
            else if ( strcmp("CM", str) == 0 ) { idx_cm = i; }
        }
    }

    if ( tr_var.fields == NULL ) return false;

    // read variant info
    plink_var_t var;
    var.schrom.assign(tr_var.str_field_at(idx_chrom));
    var.pos = tr_var.int_field_at(idx_pos);
    var.vid.assign(tr_var.str_field_at(idx_id));
    var.ref.assign(tr_var.str_field_at(idx_ref));
    split(var.alts, ",", tr_var.str_field_at(idx_alt));
    var.cM = tr_var.double_field_at(idx_cm);
    vars.push_back(var);
    ++nvars;

    return true;
}

bool PlinkReader::stream_bim() {
    if ( nvars == 0 ) {
        if ( !tr_var.open(bimf.c_str()) ) { // return false if file cannot be opened
            error("Cannot open %s", bimf.c_str());
            return false;
        } 
    }

    if ( tr_var.read_line() ) {
        plink_var_t var;
        var.schrom.assign(tr_var.str_field_at(0));
        var.vid.assign(tr_var.str_field_at(1));
        var.cM = tr_var.double_field_at(2);
        var.pos = tr_var.int_field_at(3);
        var.alts.push_back(tr_var.str_field_at(4));
        var.ref.assign(tr_var.str_field_at(5));
        vars.push_back(var);
        ++nvars;
        return true;
    }
    else {
        return false;
    }
}

bool PlinkReader::stream_bed_genos() {
    // read the variant position first
    if ( !stream_bim() ) {
        return false;
    }

    // load the genotypes if needed
    if ( !pgr_loaded ) {
        pgr.Load(bedf, (int32_t)samps.size(), samp_idx, nthreads);
        pgr_loaded = true;
    }

    // read the genotypes, read as integers
    pgr.ReadIntHardcalls(int_buf, 0, nvars-1, 0);
    return true;
}

bool PlinkReader::stream_pgen_genos() {
    // read the variant position first
    if ( !stream_pvar() ) {
        return false;
    }

    // load the genotypes if needed
    if ( !pgr_loaded ) {
        pgr.Load(pgenf, (int32_t)samps.size(), samp_idx, nthreads);
        pgr_loaded = true;
    }

    // read the genotypes, read as integers
    pgr.ReadIntHardcalls(int_buf, 0, nvars-1, 0);
    return true;
}

int32_t PlinkReader::get_variant_idx_from_cpra(const std::vector<std::string>& var_cpra, std::vector<int32_t>& variant_idx) {
    notice("bimf = %s", bimf.c_str());
    notice("pvarf = %s", pvarf.c_str());
    notice("mode_bed = %d", mode_bed);


    variant_idx.resize(var_cpra.size());
    std::fill(variant_idx.begin(), variant_idx.end(), -1);

    // create a map from CPRA to variant index
    std::map<std::string, int32_t> cpra2idx;
    for(int32_t i=0; i < var_cpra.size(); ++i) {
        cpra2idx[var_cpra[i]] = i;
    }

    // scan all variants from the 
    if ( mode_bed ) {
        return bim_cpra2idx(cpra2idx, variant_idx);
    } else {
        return pvar_cpra2idx(cpra2idx, variant_idx);
    }
}

int32_t PlinkReader::bim_cpra2idx(std::map<std::string, int32_t>& cpra2idx, std::vector<int32_t>& variant_idx) {
    tsv_reader tr_bim(bimf.c_str());
    char buf[65535];
    int32_t nv = 0, nmatch = 0;
    while ( tr_bim.read_line() ) {
        snprintf(buf, 65535, "%s:%s:%s:%s", tr_bim.str_field_at(0), tr_bim.str_field_at(3), tr_bim.str_field_at(5), tr_bim.str_field_at(4));
        if ( tr_bim.int_field_at(3) == 271129 ) 
            notice("buf = %s, cpra2idx.size() = %zu", buf, cpra2idx.size());
        std::map<std::string, int32_t>::iterator it = cpra2idx.find(buf);
        if ( it != cpra2idx.end() ) {
            variant_idx[it->second] = nv;
            ++nmatch;
        }
        ++nv;
    }
    tr_bim.close();
    return nmatch;
}

int32_t PlinkReader::pvar_cpra2idx(std::map<std::string, int32_t>& cpra2idx, std::vector<int32_t>& variant_idx) {
    tsv_reader tr_pvar(pvarf.c_str());
    char buf[65535];
    int32_t nv = 0, nmatch = 0;
    int32_t idx_chrom = -1, idx_pos = -1, idx_id = -1, idx_ref = -1, idx_alt = -1, idx_qual = -1, idx_cm = -1, idx_filter = -1, idx_format = -1, idx_info = -1;
    while ( tr_pvar.read_line() ) {
        const char* s = tr_pvar.str_field_at(0);
        if ( s[0] == '#' ) {
            if ( strncmp(s, "#CHROM", 6) == 0 ) {  // meta line
                // this is a header line that contains column names
                for(int32_t i=0; i < tr_pvar.nfields; ++i) {
                    const char* str = tr_pvar.str_field_at(i) + ( i == 0 ? 1 : 0 );
                    if ( strcmp("CHROM", str) == 0 ) { idx_chrom = i; }
                    else if ( strcmp("POS", str) == 0 ) { idx_pos = i; }
                    else if ( strcmp("ID", str) == 0 ) { idx_id = i; }
                    else if ( strcmp("REF", str) == 0 ) { idx_ref = i; }
                    else if ( strcmp("ALT", str) == 0 ) { idx_alt = i; }
                    else if ( strcmp("QUAL", str) == 0 ) { idx_qual = i; }
                    else if ( strcmp("FILTER", str) == 0 ) { idx_filter = i; }
                    else if ( strcmp("INFO", str) == 0 ) { idx_info = i; }
                    else if ( strcmp("FORMAT", str) == 0 ) { idx_format = i; }
                    else if ( strcmp("CM", str) == 0 ) { idx_cm = i; }
                }
                if ( ( idx_chrom < 0 ) || ( idx_pos < 0 ) || ( idx_ref < 0 ) || ( idx_alt < 0 ) ) {
                    error("Cannot find the required columns in %s", pvarf.c_str());
                }
            }
        }
        else {
            snprintf(buf, 65535, "%s:%s:%s:%s", tr_pvar.str_field_at(idx_chrom), tr_pvar.str_field_at(idx_pos), tr_pvar.str_field_at(idx_ref), tr_pvar.str_field_at(idx_alt));
            //error("%s %d %d %d %d", buf, idx_chrom, idx_pos, idx_ref, idx_alt);
            std::map<std::string, int32_t>::iterator it = cpra2idx.find(buf);
            if ( it != cpra2idx.end() ) {
                variant_idx[it->second] = nv;
                ++nmatch;
            }
            ++nv;
        }
    }
    tr_pvar.close();
    return nmatch;
}

bool PlinkReader::get_bed_genos_at(int32_t var_idx) {
    // load the genotypes if needed
    if ( !pgr_loaded ) {
        pgr.Load(bedf, (int32_t)samps.size(), samp_idx, nthreads);
        pgr_loaded = true;
    }

    // read the genotypes, read as integers
    pgr.ReadIntHardcalls(int_buf, 0, var_idx, 0);
    return true;
}

bool PlinkReader::get_pgen_genos_at(int32_t var_idx) {
    // load the genotypes if needed
    if ( !pgr_loaded ) {
        pgr.Load(pgenf, (int32_t)samps.size(), samp_idx, nthreads);
        pgr_loaded = true;
    }

    // read the genotypes, read as integers
    pgr.ReadIntHardcalls(int_buf, 0, var_idx, 0);
    return true;
}
