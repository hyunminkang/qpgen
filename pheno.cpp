#include "pheno.h"
#include "qgenlib/qgen_error.h"
#include "qgenlib/tsv_reader.h"
#include "qpgen_utils.h"

#include <cstring>

int32_t PhenoMatrix::subset_sample_ids(const std::vector<std::string>& samp_ids) {
    // identify overlapping indices
    std::vector<int32_t> phe_sample_indices;
    std::vector<int32_t> given_sample_indices;
    int32_t n_overlaps = index_overlapping_ids(this->samp_ids, samp_ids, phe_sample_indices, given_sample_indices);
    if ( n_overlaps == 0 ) {
        error("No overlapping sample IDs found between phenotype matrix and given sample IDs");
    }
    // subset the phenotype matrix
    Eigen::MatrixXd pheno_mat_sub(n_overlaps, pheno_mat.cols());
    Eigen::Vector<bool, Eigen::Dynamic> pheno_mask_sub(n_overlaps, pheno_mat.cols());
    std::vector<std::string> samp_ids_sub(n_overlaps);
    for ( int32_t i = 0; i < n_overlaps; ++i ) {
        pheno_mat_sub.row(i) = pheno_mat.row( phe_sample_indices[i] );
        pheno_mask_sub.row(i) = pheno_mask.row( phe_sample_indices[i] );
        samp_ids_sub[i] = this->samp_ids[ phe_sample_indices[i] ];
    }
    pheno_mat = pheno_mat_sub;
    pheno_mask = pheno_mask_sub;
    this->samp_ids = samp_ids_sub;
    return n_overlaps;
} 

int32_t PhenoMatrix::subset_pheno_ids(const std::vector<std::string>& pheno_ids) {
    // identify overlapping indices
    std::vector<int32_t> phe_pheno_indices;
    std::vector<int32_t> given_pheno_indices;
    int32_t n_overlaps = index_overlapping_ids(this->pheno_ids, pheno_ids, phe_pheno_indices, given_pheno_indices);
    //notice("Found %d overlapping phenotype IDs between phenotype matrix and given phenotype IDs", n_overlaps);
    if ( n_overlaps == 0 ) {
        error("No overlapping phenotype IDs found between phenotype matrix and given phenotype IDs");
    }
    // subset the phenotype matrix
    Eigen::MatrixXd pheno_mat_sub(pheno_mat.rows(), n_overlaps);
    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> pheno_mask_sub(pheno_mat.rows(), n_overlaps);
    std::vector<std::string> pheno_ids_sub(n_overlaps);
    std::vector<genomeLocus> pheno_loci_sub;
    for ( int32_t i = 0; i < n_overlaps; ++i ) {
        pheno_mat_sub.col(i) = pheno_mat.col( phe_pheno_indices[i] );
        pheno_mask_sub.col(i) = pheno_mask.col( phe_pheno_indices[i] );
        pheno_ids_sub[i] = this->pheno_ids[ phe_pheno_indices[i] ];
        if ( has_loci ) {
            pheno_loci_sub.push_back( this->pheno_loci[ phe_pheno_indices[i] ] );
        }
    }
    pheno_mat = pheno_mat_sub;
    pheno_mask = pheno_mask_sub;
    this->pheno_ids = pheno_ids_sub;
    if ( has_loci ) {
        this->pheno_loci = pheno_loci_sub;
        this->pheno_locusmap.clear();
        for ( int32_t i = 0; i < (int32_t)pheno_loci_sub.size(); ++i ) {
            const genomeLocus& locus = pheno_loci_sub[i];
            this->pheno_locusmap.add( locus.chrom.c_str(), locus.beg1, locus.end0, i );
        }
    }
    return n_overlaps;
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
        return true;
    }
    else {
        error("Unsupported phenotype format: %s", pheno_format);
        return false;
    }
}
