#ifndef TAGSORT_PARTIAL_FILE_MERGE_H_
#define TAGSORT_PARTIAL_FILE_MERGE_H_

#include <string>
#include <vector>

#include "tagsort_input_options.h"

// Takes the sorted partial files produced by splitAndPartialSortToFiles(),
// and merges them into a single sorted output file.
// Returns the number of alignments written into that output file.
int mergePartialFiles(INPUT_OPTIONS_TAGSORT const& options,
                      std::vector<std::string> const& partial_files);

#endif // TAGSORT_PARTIAL_FILE_MERGE_H_
