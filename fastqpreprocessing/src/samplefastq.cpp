#include "fastq_common.h"
#include "input_options.h"
#include <fstream>

std::vector<std::pair<char, int>> parseReadStructure(std::string const& read_structure)
{
  std::vector<std::pair<char, int>> ret;
  int next_ind = 0;
  while (next_ind < read_structure.size())
  {
    int type_ind = read_structure.find_first_not_of("0123456789", next_ind);
    assert(type_ind != std::string::npos);
    char type = read_structure[type_ind];
    int len = std::stoi(read_structure.substr(next_ind, type_ind - next_ind));
    ret.emplace_back(type, len);
    next_ind = type_ind + 1;
  }
  return ret;
}

std::vector<std::pair<char, int>> g_parsed_read_structure;

void fillSamRecordWithReadStructure(SamRecord* sam, FastQFile* fastQFileI1,
                                    FastQFile* fastQFileR1, FastQFile* fastQFileR2, FastQFile* fastQFileR3,
                                    bool has_I1_file_list, bool has_R3_file_list,
                                    std::string orientation)
{
  // check the sequence names matching
  std::string a = std::string(fastQFileR1->myRawSequence.c_str());
  std::string b = std::string(fastQFileR1->myQualityString.c_str());
  
  // extract the raw barcode and UMI 8C18X6C9M1X and raw barcode and UMI quality string
  std::string barcode_seq, barcode_quality, umi_seq, umi_quality;
  int cur_ind = 0;
  
  for (auto [tag, length] : g_parsed_read_structure)
  {
    switch (tag)
    {
    case 'C':
      barcode_seq += a.substr(cur_ind, length);
      barcode_quality += b.substr(cur_ind, length);
      break;
    case 'M':
      umi_seq += a.substr(cur_ind, length);
      umi_quality += b.substr(cur_ind, length);
      break;
    default:
      break;
    }
    cur_ind += length;
  }
  fillSamRecordCommon(sam, fastQFileI1, fastQFileR1, fastQFileR2, 
                      fastQFileR3, has_I1_file_list, has_R3_file_list,
                      barcode_seq, barcode_quality, umi_seq, umi_quality); 
}

std::string slideseqBarcodeGetter(SamRecord* sam)
{
  return std::string(sam->getString("CR").c_str());
}

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

  g_parsed_read_structure = parseReadStructure(options.read_structure);
  mainCommon(options.white_list_file, options.barcode_orientation, /*num_writer_threads=*/1, options.output_format,
             options.I1s, options.R1s, options.R2s, options.R3s, options.sample_id,
             fillSamRecordWithReadStructure, slideseqBarcodeGetter,
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
