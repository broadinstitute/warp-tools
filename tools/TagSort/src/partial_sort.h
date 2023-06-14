#ifndef TAGSORT_PARTIAL_SORT_H_
#define TAGSORT_PARTIAL_SORT_H_

#include "tagsort_input_options.h"

#include <string>
#include <vector>

// From the input bam, arbitrarily split the records into a number of groups,
// sort each group, and write each sorted group to a file, each randomly named.
// Returns the file paths of the created sorted partial files.
std::vector<std::string> splitAndPartialSortToFiles(INPUT_OPTIONS_TAGSORT options);

#endif // TAGSORT_PARTIAL_SORT_H_
