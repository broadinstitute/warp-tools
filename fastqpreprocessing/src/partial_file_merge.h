#ifndef TAGSORT_PARTIAL_FILE_MERGE_H_
#define TAGSORT_PARTIAL_FILE_MERGE_H_

#include <string>
#include <vector>

#include "tagsort_input_options.h"

// returns number of alignments processed
int mergePartialFiles(INPUT_OPTIONS_TAGSORT const& options,
                      std::vector<std::string> const& partial_files);

#endif // TAGSORT_PARTIAL_FILE_MERGE_H_
