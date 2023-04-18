#include "fastq_common.h"
#include "input_options.h"

void outputHandler(WriteQueue* cur_write_queue, SamRecord* samrec, int reader_thread_index)
{
  cur_write_queue->enqueueWrite(std::make_pair(samrec, reader_thread_index));
}

int main(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options = readOptionsFastqSlideseq(argc, argv);
  // number of output bam files, and one writer thread per bam file
  int num_writer_threads = get_num_blocks(options);

  std::vector<std::pair<char, int>> g_parsed_read_structure;
  g_parsed_read_structure = parseReadStructure(options.read_structure);

  mainCommon(options.white_list_file, options.barcode_orientation, num_writer_threads, options.output_format,
             options.I1s, options.R1s, options.R2s, options.R3s, options.sample_id, g_parsed_read_structure,
             //fillSamRecordWithReadStructure, slideseqBarcodeGetter, 
             outputHandler);

  return 0;
}
