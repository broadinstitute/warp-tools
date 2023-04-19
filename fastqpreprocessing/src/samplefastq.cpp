#include "fastq_common.h"
#include "input_options.h"
#include <fstream>

void outputHandler(WriteQueue* cur_write_queue, SamRecord* samrec, int reader_thread_index)
{
  cur_write_queue->enqueueWrite(std::make_pair(samrec, reader_thread_index));
}

int main(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options = readOptionsFastqSlideseq(argc, argv);

  std::ofstream outfile_r1("sampled_down.R1");
  if (!outfile_r1)
    crash("Failed to open output file sampled_down.R1");
  std::ofstream outfile_r2("sampled_down.R2");
  if (!outfile_r2)
    crash("Failed to open output file sampled_down.R2");

  std::vector<std::pair<char, int>> g_parsed_read_structure;
  g_parsed_read_structure = parseReadStructure(options.read_structure);

  mainCommon(options.white_list_file, options.barcode_orientation, /*num_writer_threads=*/1, options.output_format,
             options.I1s, options.R1s, options.R2s, options.R3s, options.sample_id, g_parsed_read_structure,
             //fillSamRecordWithReadStructure, slideseqBarcodeGetter, 
             [&outfile_r1, &outfile_r2](WriteQueue* ignored1, SamRecord* sam, int reader_thread_index)
             {
               if (sam->getStringTag("CB"))
               {
                   // Read structure can be 8C18X6C9M1X, 16C10M or 16C12M 
                   const char* barcode = sam->getString("CR").c_str();
                   const char* quality_score = sam->getString("CY").c_str();
                  
                   // Barcode 
                   int barcode_index = 0;        
                   outfile_r1 << "@" <<  sam->getReadName() << "\n";
                   for (auto [tag, length] : g_parsed_read_structure)
                   {
                      if (tag == 'C')
                      {
                            outfile_r1 << std::string_view(barcode + barcode_index, length);
                            barcode_index = length;
                      }
                      else if (tag == 'M')
                            outfile_r1 << sam->getString("UR");
                      else if (tag == 'X')
                            outfile_r1 << std::string_view("CTTCAGCGTTCCCGAGAG", length);        
                   } 
                   outfile_r1 << "\n" << "+\n";
                   
                   // Quality 
                   barcode_index = 0;
                   for (auto [tag, length] : g_parsed_read_structure)
                   {
                      if (tag == 'C')
                      {
                            outfile_r1 << std::string_view(quality_score + barcode_index, length);
                            barcode_index = length;
                      }
                      else if (tag == 'M')
                            outfile_r1 << sam->getString("UY");
                      else if (tag == 'X')
                            outfile_r1 << std::string_view("FFFFFFFFFFFFFFFFFF", length);
                   }       
                   outfile_r1 << "\n";
                 
                   outfile_r2 << "@" << sam->getReadName() << "\n"
                          << sam->getSequence() << "\n"
                          << "+\n"
                          << sam->getQuality() << "\n";
               }
               releaseReaderThreadMemory(reader_thread_index,sam);
             });
  return 0;
}
