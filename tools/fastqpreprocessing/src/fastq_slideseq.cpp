#include "fastq_common.h"
#include "input_options.h"

int main(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options = readOptionsFastqSlideseq(argc, argv);
  // number of output bam files, and one writer thread per bam file
  int num_writer_threads = get_num_blocks(options);
  // hardcoded this to 1000 in case of large files
  num_writer_threads =  (num_writer_threads > 1000) ? 1000 : num_writer_threads;

  std::vector<std::pair<char, int>> g_parsed_read_structure = parseReadStructure(options.read_structure);

  mainCommon(options.white_list_file, options.barcode_orientation, num_writer_threads, options.output_format,
             options.I1s, options.R1s, options.R2s, options.R3s, options.sample_id, g_parsed_read_structure,
             options.sample_bool);

  return 0;
}
