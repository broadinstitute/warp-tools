#include "fastq_common.h"
#include "input_options.h"

int main(int argc, char** argv)
{
  InputOptionsFastqProcess options = readOptionsFastqProcess(argc, argv);
  // number of output bam files, and one writer thread per bam file
  int num_writer_threads = get_num_blocks(options);
  // hardcoded this to 1000 in case of large files
  num_writer_threads =  (num_writer_threads > 1000) ? 1000 : num_writer_threads;

  // added this for consistency with other code
  std::vector<std::pair<char, int>> g_parsed_read_structure = parseReadStructure(options.read_structure);

  mainCommon(options.white_list_file, options.barcode_orientation, num_writer_threads, options.output_format,
             options.I1s, options.R1s, options.R2s, options.R3s, options.sample_id, g_parsed_read_structure,
             options.sample_bool);

  return 0;
}

int main(int argc, char** argv)
{
  InputOptionsFastqProcess options = readOptionsFastqProcess(argc, argv);

  if (options.bam_size == -1.0) { // Check if bam_size is still at its uninitialized state
    std::vector<std::string> allFiles;
    allFiles.insert(allFiles.end(), options.I1s.begin(), options.I1s.end());
    allFiles.insert(allFiles.end(), options.R1s.begin(), options.R1s.end());
    allFiles.insert(allFiles.end(), options.R2s.begin(), options.R2s.end());
    allFiles.insert(allFiles.end(), options.R3s.begin(), options.R3s.end());

    int64_t totalSize = calculateTotalSize(allFiles);

    // Convert total size to GB if it is in bytes
    options.bam_size = static_cast<double>(totalSize) / (1024 * 1024 * 1024);
  }

  // number of output bam files, and one writer thread per bam file
  int num_writer_threads = get_num_blocks(options);

  // hardcoded this to 1000 in case of large files
  num_writer_threads = (num_writer_threads > 1000) ? 1000 : num_writer_threads;

  // added this for consistency with other code
  std::vector<std::pair<char, int>> g_parsed_read_structure = parseReadStructure(options.read_structure);

  mainCommon(options.white_list_file, options.barcode_orientation, num_writer_threads, options.output_format,
             options.I1s, options.R1s, options.R2s, options.R3s, options.sample_id, g_parsed_read_structure,
             options.sample_bool);

  return 0;
}
