#include "pheno.h"
#include "qgenlib/qgen_error.h"
#include "qgenlib/tsv_reader.h"
#include "qpgen_utils.h"

#include <cstring>

bool PhenoMatrix::subset_sample_pheno_indices(const std::vector<int32_t>& samp_indices, const std::vector<int32_t>& pheno_indices) {
    if ( samp_indices.size() == this->samp_ids.size() && pheno_indices.size() == this->pheno_ids.size() ) {
        bool change_needed = false;
        for ( int i = 0; i < samp_indices.size(); ++i ) {
            if ( samp_indices[i] != i ) {
                change_needed = true;
                break;
            }
        }
        if ( !change_needed ) {
            for ( int j = 0; j < pheno_indices.size(); ++j ) {
                if ( pheno_indices[j] != j ) {
                    change_needed = true;
                    break;
                }
            }
        }
        if ( !change_needed ) {
            return true;
        }
    }

    std::vector<std::string> samp_ids_sub;
    std::vector<std::string> pheno_ids_sub;
    Eigen::MatrixXd pheno_mat_sub(samp_indices.size(), pheno_indices.size());
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> pheno_mask_sub(samp_indices.size(), pheno_indices.size());
    std::vector<genomeLocus> pheno_loci_sub;
    for ( size_t i = 0; i < samp_indices.size(); ++i ) {
        int32_t sidx = samp_indices[i];
        if ( sidx < 0 || sidx >= (int32_t)samp_ids.size() ) {
            error("Sample index %d is out of bounds [0,%d)", sidx, (int32_t)samp_ids.size());
        }
        samp_ids_sub.push_back( samp_ids[sidx] );
        for ( size_t j = 0; j < pheno_indices.size(); ++j ) {
            int32_t pidx = pheno_indices[j];
            if ( pidx < 0 || pidx >= (int32_t)this->pheno_ids.size() ) {
                error("Phenotype index %d is out of bounds [0,%d)", pidx, (int32_t)this->pheno_ids.size());
            }
            if ( i == 0 ) {
                pheno_ids_sub.push_back( this->pheno_ids[pidx] );
                if ( has_loci ) {
                    pheno_loci_sub.push_back( this->pheno_loci[pidx] );
                }
            }
            pheno_mat_sub(i,j) = pheno_mat(sidx, pidx);
            pheno_mask_sub(i,j) = pheno_mask(sidx, pidx);
        }
    }
    pheno_mat = pheno_mat_sub;
    pheno_mask = pheno_mask_sub;
    this->samp_ids = samp_ids_sub;
    this->pheno_ids = pheno_ids_sub;
    if ( has_loci ) {
        this->pheno_loci = pheno_loci_sub; 
    }
    rebuild_id2index_map(samp_ids, samp_id2idx);
    rebuild_id2index_map(pheno_ids, pheno_id2idx);
    return true;
}

int32_t PhenoMatrix::subset_pheno_indices(const std::vector<int32_t>& pheno_indices) {
    if ( pheno_indices.size() == this->pheno_ids.size() ) {
        bool change_needed = false;
        for ( int j = 0; j < pheno_indices.size(); ++j ) {
            if ( pheno_indices[j] != j ) {
                change_needed = true;
                break;
            }
        }
        if ( !change_needed ) {
            return (int32_t)this->pheno_ids.size();
        }
    }

    std::vector<std::string> pheno_ids_sub;
    Eigen::MatrixXd pheno_mat_sub(pheno_mat.rows(), pheno_indices.size());
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> pheno_mask_sub(pheno_mat.rows(), pheno_indices.size());
    std::vector<genomeLocus> pheno_loci_sub;
    for ( size_t j = 0; j < pheno_indices.size(); ++j ) {
        int32_t idx = pheno_indices[j];
        if ( idx < 0 || idx >= (int32_t)pheno_ids.size() ) {
            error("Phenotype index %d is out of bounds [0,%d)", idx, (int32_t)this->pheno_ids.size());
        }
        pheno_ids_sub.push_back( pheno_ids[idx] );
        if ( has_loci ) {
            pheno_loci_sub.push_back( pheno_loci[idx] );
        }
        pheno_mat_sub.col(j) = pheno_mat.col( idx );
        pheno_mask_sub.col(j) = pheno_mask.col( idx );
    }
    pheno_mat = pheno_mat_sub;
    pheno_mask = pheno_mask_sub;
    this->pheno_ids = pheno_ids_sub;
    if ( has_loci ) {
        this->pheno_loci = pheno_loci_sub;
    }

    rebuild_id2index_map(pheno_ids, pheno_id2idx);

    return (int32_t)this->pheno_ids.size();
}

int32_t PhenoMatrix::subset_sample_indices(const std::vector<int32_t>& samp_indices) {
    if ( samp_indices.size() == this->samp_ids.size() ) {
        bool change_needed = false;
        for ( int i = 0; i < samp_indices.size(); ++i ) {
            if ( samp_indices[i] != i ) {
                change_needed = true;
                break;
            }
        }
        if ( !change_needed ) {
            return (int32_t)this->samp_ids.size();
        }
    }

    std::vector<std::string> samp_ids_sub;
    Eigen::MatrixXd pheno_mat_sub(samp_indices.size(), pheno_mat.cols());
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> pheno_mask_sub(samp_indices.size(), pheno_mat.cols());
    for ( size_t i = 0; i < samp_indices.size(); ++i ) {
        int32_t idx = samp_indices[i];
        if ( idx < 0 || idx >= (int32_t)samp_ids.size() ) {
            error("Sample index %d is out of bounds [0,%d)", idx, (int32_t)samp_ids.size());
        }
        samp_ids_sub.push_back( samp_ids[idx] );
        pheno_mat_sub.row(i) = pheno_mat.row( idx );
        pheno_mask_sub.row(i) = pheno_mask.row( idx );
    }
    pheno_mat = pheno_mat_sub;
    pheno_mask = pheno_mask_sub;
    this->samp_ids = samp_ids_sub;

    rebuild_id2index_map(samp_ids, samp_id2idx);

    return (int32_t)this->samp_ids.size();
}

bool PhenoMatrix::subset_sample_pheno_ids(const std::vector<std::string>& samp_ids, const std::vector<std::string>& pheno_ids) {
    std::vector<int32_t> phe_sample_indices;
    std::vector<int32_t> phe_pheno_indices;
    for(int32_t i = 0; i < samp_ids.size(); ++i ) {
        const std::string& sid = samp_ids[i];
        std::map<std::string, int32_t>::const_iterator it = samp_id2idx.find( sid );
        if ( it != samp_id2idx.end() ) {
            phe_sample_indices.push_back( it->second );
        }
        else {
            error("Sample ID '%s' not found in phenotype matrix", sid.c_str());
        }
    }
    for(int32_t i = 0; i < pheno_ids.size(); ++i ) {
        const std::string& pid = pheno_ids[i];
        std::map<std::string, int32_t>::const_iterator it = pheno_id2idx.find( pid );
        if ( it != pheno_id2idx.end() ) {
            phe_pheno_indices.push_back( it->second );
        }
        else {
            error("Phenotype ID '%s' not found in phenotype matrix", pid.c_str());
        }
    }
    return subset_sample_pheno_indices(phe_sample_indices, phe_pheno_indices);
}


int32_t PhenoMatrix::subset_sample_ids(const std::vector<std::string>& samp_ids) {
    std::vector<int32_t> phe_sample_indices;
    for(int32_t i = 0; i < samp_ids.size(); ++i ) {
        const std::string& sid = samp_ids[i];
        std::map<std::string, int32_t>::const_iterator it = samp_id2idx.find( sid );
        if ( it != samp_id2idx.end() ) {
            phe_sample_indices.push_back( it->second );
        }
        else {
            error("Sample ID '%s' not found in phenotype matrix", sid.c_str());
        }
    }
    return subset_sample_indices( phe_sample_indices );
} 

int32_t PhenoMatrix::subset_pheno_ids(const std::vector<std::string>& pheno_ids) {
    std::vector<int32_t> phe_pheno_indices;
    for(int32_t i = 0; i < pheno_ids.size(); ++i ) {
        const std::string& pid = pheno_ids[i];
        std::map<std::string, int32_t>::const_iterator it = pheno_id2idx.find( pid );
        if ( it != pheno_id2idx.end() ) {
            phe_pheno_indices.push_back( it->second );
        }
        else {
            error("Phenotype ID '%s' not found in phenotype matrix", pid.c_str());
        }
    }
    return subset_pheno_indices( phe_pheno_indices );
}

bool PhenoMatrix::load_pheno_matrix(const char* pheno_file, const char* pheno_format, const char delim) {
    format.assign(pheno_format);
    if ( format.compare("regenie") == 0 ||
         format.compare("tsv-sample-row") == 0 ) { // sample in rows, phenos in columns
        notice("Loading sample-row format phenotype matrix from %s", pheno_file);
        int32_t icol_indid = -1;
        int32_t icol_data = -1;
        if ( format.compare("regenie") == 0 ) {
            icol_indid = 1;
            icol_data = 2;
        }
        else {
            icol_indid = 0;
            icol_data = 1;
        }
        tsv_reader tr(pheno_file);
        if ( delim != '\0' ) {
            tr.delimiter = delim; // note that defeault delimiter is whitespace
        }
        // read header
        int32_t nlines = 0;
        int32_t nfields = -1;
        int32_t nrow_est = 100;
        int32_t ntraits = -1;
        int32_t nsamps = 0;
        while ( tr.read_line() ) {
            if ( nlines == 0 ) {
                nfields = tr.nfields;
                for ( int32_t i = icol_data; i < nfields; ++i ) {
                    pheno_ids.push_back( std::string(tr.str_field_at(i)) );
                }
                ntraits = nfields - icol_data;
                notice("Loading %d traits from phenotype matrix", ntraits);
                pheno_mat.resize(nrow_est, ntraits);
                pheno_mask.resize(nrow_est, ntraits);
                ++nlines;
            }
            else {
                if ( nfields != tr.nfields ) {
                    error("Inconsistent number of fields at line %d in phenotype file %s", nlines + 1, pheno_file);
                }
                // check if expansion is needed
                if ( nsamps >= nrow_est ) {
                    nrow_est = nrow_est * 2;
                    notice("Expanding phenotype matrix to %d samples", nrow_est);
                    pheno_mat.conservativeResize(nrow_est, ntraits);
                    pheno_mask.conservativeResize(nrow_est, ntraits);
                }
                samp_ids.push_back( std::string(tr.str_field_at(icol_indid)) );
                for ( int32_t i = icol_data; i < nfields; ++i ) {
                    if ( is_missing( tr.str_field_at(i) ) ) {
                        pheno_mask(nsamps, i - icol_data) = false;
                        pheno_mat(nsamps, i - icol_data) = 0;
                        has_missing = true;
                    }
                    else {
                        pheno_mask(nsamps, i - icol_data) = true;
                        pheno_mat(nsamps, i - icol_data) = tr.double_field_at(i);
                    }
                }
                ++nlines;
                ++nsamps;
            }
        }
        notice("Resizing phenotype matrix for %d samples and %d traits", nsamps, ntraits);
        // resize to the actual number of samples
        pheno_mat.conservativeResize(nsamps, ntraits);
        pheno_mask.conservativeResize(nsamps, ntraits);

        rebuild_id2index_map(samp_ids, samp_id2idx);
        rebuild_id2index_map(pheno_ids, pheno_id2idx);
        return true;
    }
    else if ( format.compare("tsv-sample-col") == 0 ||
              format.compare("tensorqtl") == 0 ) {
        // sample in columns, phenos in rows
        int32_t icol_pheid = -1;
        int32_t icol_chrom = -1;
        int32_t icol_beg = -1;
        int32_t icol_end = -1;
        int32_t icol_data = -1;
        if ( format.compare("tensorqtl") == 0 ) {
            icol_chrom = 0;
            icol_beg = 1;
            icol_end = 2;
            has_loci = true;
            icol_pheid = 3;
            icol_data = 4;
        }
        else {
            icol_pheid = 0;
            icol_data = 1;
        }
        tsv_reader tr(pheno_file);
        if ( delim != '\0' ) {
            tr.delimiter = delim; // note that defeault delimiter is whitespace
        }
        // read header
        int32_t nlines = 0;
        int32_t nfields = -1;
        int32_t ncol_est = 100;
        int32_t nsamps = -1;
        int32_t ntraits = 0;
        while ( tr.read_line() ) {
            if ( nlines == 0 ) {
                nfields = tr.nfields;
                if ( has_loci ) {
                    if ( strcmp("#chr", tr.str_field_at(0)) != 0 ||
                         strcmp("start", tr.str_field_at(1)) != 0 ||
                         strcmp("end", tr.str_field_at(2)) != 0 ) {
                        error("Expected chr,start,end columns for TensorQTL format in pheno file %s", pheno_file);
                    }
                }
                for ( int32_t i = icol_data; i < nfields; ++i ) {
                    samp_ids.push_back( std::string(tr.str_field_at(i)) );
                }
                nsamps = nfields - icol_data;
                pheno_mat.resize(nsamps, ncol_est);
                pheno_mask.resize(nsamps, ncol_est);
                nlines++;
            }
            else {
                if ( nfields != tr.nfields ) {
                    error("Inconsistent number of fields at line %d in phenotype file %s", nlines + 1, pheno_file);
                }
                // check if expansion is needed
                if ( ntraits >= ncol_est ) {
                    ncol_est = ncol_est * 2;
                    pheno_mat.conservativeResize(nsamps, ncol_est);
                    pheno_mask.conservativeResize(nsamps, ncol_est);
                }
                pheno_ids.push_back( std::string(tr.str_field_at(icol_pheid)) );
                if ( has_loci ) {
                    const char* chrom = tr.str_field_at(icol_chrom);
                    int32_t beg1 = tr.int_field_at(icol_beg);
                    int32_t end0 = tr.int_field_at(icol_end);
                    pheno_loci.push_back( genomeLocus(chrom, beg1, end0) );
                    pheno_locusmap.add(chrom, beg1, end0, ntraits);
                }
                for ( int32_t i = icol_data; i < nfields; ++i ) {
                    if ( is_missing( tr.str_field_at(i) ) ) {
                        pheno_mask(i - icol_data, ntraits) = false;
                        pheno_mat(i - icol_data, ntraits) = 0;
                        has_missing = true;
                    }
                    else {
                        pheno_mask(i - icol_data, ntraits) = true;
                        pheno_mat(i - icol_data, ntraits) = tr.double_field_at(i);
                    }
                }
                ++nlines;
                ++ntraits;
            }
        }
        // resize to the actual number of samples
        pheno_mat.conservativeResize(nsamps, ntraits);
        pheno_mask.conservativeResize(nsamps, ntraits);

        rebuild_id2index_map(samp_ids, samp_id2idx);
        rebuild_id2index_map(pheno_ids, pheno_id2idx);
        return true;
    }
    else {
        error("Unsupported phenotype format: %s", pheno_format);
        return false;
    }
}

bool PhenoMatrix::sample_ids_sorted() const {
    for ( size_t i = 1; i < samp_ids.size(); ++i ) {
        if ( samp_ids[i] < samp_ids[i-1] ) {
            return false;
        }
    }
    return true;
}

bool PhenoMatrix::pheno_ids_sorted() const {
    for ( size_t i = 1; i < pheno_ids.size(); ++i ) {
        if ( pheno_ids[i] < pheno_ids[i-1] ) {
            return false;
        }
    }
    return true;
}