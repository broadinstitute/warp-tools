#include "fastq_common.h"
#include "input_options.h"

int main(int argc, char** argv)
{
  InputOptionsFastqProcess options = readOptionsFastqProcess(argc, argv);
  
  int num_writer_threads = 1;
  if (options.num_output_files != -1) {
    std::cout<<"options.num_output_file != -1\n";
    num_writer_threads = options.num_output_files;
  }
  else {
    std::cout<<"ELSE\n";
    // number of output bam files, and one writer thread per bam file
    num_writer_threads = get_num_blocks(options);
    // hardcoded this to 1000 in case of large files
    num_writer_threads =  (num_writer_threads > 1000) ? 1000 : num_writer_threads;
  }

  std::cout<<num_writer_threads<<"\n";
  // num_writer_threads= (options.num_output_files != -1) ? options.num_output_files : num_writer_threads;
  std::cout<<options.num_output_files<<"\n";

  // added this for consistency with other code
  std::vector<std::pair<char, int>> g_parsed_read_structure = parseReadStructure(options.read_structure);

  mainCommon(options.white_list_file, options.barcode_orientation, num_writer_threads, options.output_format,
             options.I1s, options.R1s, options.R2s, options.R3s, options.sample_id, g_parsed_read_structure,
             options.sample_bool);

  return 0;
}
