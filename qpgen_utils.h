#ifndef __QPGEN_UTILS_H
#define __QPGEN_UTILS_H

#include <vector>
#include <string>

#include "Eigen/Dense"
#include "qgenlib/qgen_error.h"
#include "qgenlib/hts_utils.h"
#include "qpgen.h"

// identify overlapping IDS from list1 and list2, and store sorted overlapping IDs in overlapping_ids
int32_t identify_overlapping_ids(const std::vector<std::string>& list1, const std::vector<std::string>& list2, std::vector<std::string>& overlapping_ids);

// index overlapping IDs from list1 and list2, and store sorted indices in index1 and index2
int32_t index_overlapping_ids(const std::vector<std::string>& list1, const std::vector<std::string>& list2, std::vector<int32_t>& index1, std::vector<int32_t>& index2);

void fill_list_indices(const std::vector<std::string>& list, const std::map<std::string, int32_t>& list2idx, std::vector<int32_t>& indices);

int32_t rebuild_id2index_map(
    const std::vector<std::string>& ids,
    std::map<std::string, int32_t>& id2idx_map
);

#endif // __QPGEN_UTILS_H