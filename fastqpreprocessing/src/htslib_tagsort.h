#ifndef __HTSLIB_TAG_SORT__
#define __HTSLIB_TAG_SORT__

#include <htslib/sam.h>
#include "input_options.h"
#include "utilities.h"


// From the input bam, arbitrarily split the records into a number of groups,
// sort each group, and write each sorted group to a file, each randomly named.
// Returns the file paths of the created sorted partial files.
std::vector<std::string> create_sorted_file_splits_htslib(INPUT_OPTIONS_TAGSORT options);

#endif
