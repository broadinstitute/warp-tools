#include "fastq_common.h"
#include "input_options.h"

// ---------------------------------------------------
// Global Variables
// ----------------------------------------------------
std::vector<std::pair<char, int>> g_parsed_read_structure;

// ---------------------------------------------------
// Main
// ----------------------------------------------------
int main(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options = readOptionsFastqSlideseq(argc, argv);
  
  g_parsed_read_structure = parseReadStructure(options.read_structure);

  mainCommon(options.white_list_file, options.barcode_orientation, /*num_writer_threads=*/1, options.output_format,
             options.I1s, options.R1s, options.R2s, options.R3s, options.sample_id, g_parsed_read_structure,
             options.sample_bool);
  return 0;
}
