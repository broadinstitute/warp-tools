#ifndef __SCTOOLS_FASTQPREPROCESSING_INPUT_OPTIONS_H_
#define __SCTOOLS_FASTQPREPROCESSING_INPUT_OPTIONS_H_

#include <string>
#include <vector>

void crash(std::string msg);

struct INPUT_OPTIONS_FASTQ_READ_STRUCTURE
{
  // I1, R1 and R2 files name
  std::vector<std::string> I1s, R1s, R2s, R3s;

  // Bead Barcode list
  std::string white_list_file;

   // Barcode orientation 
  std::string barcode_orientation = "FIRST_BP";

  std::string output_format;

  // Bam file size to split by (in GB)
  double bam_size = 1.0;

  std::string read_structure;

  std::string sample_id;

  // if set to false we print out all valid/invalid barcodes.
  bool sample_bool = false;
};

// Structure to hold input options for fastqprocess
struct InputOptionsFastqProcess
{
  // I1, R1 and R2 files name
  std::vector<std::string> I1s, R1s, R2s, R3s;

  // Barcode white list file
  std::string white_list_file;

  // Barcode orientation 
  std::string barcode_orientation = "FIRST_BP";

  std::string output_format;

  // Bam file size to split by (in GB)
  double bam_size = 1.0;

  std::string read_structure;

  std::string sample_id;

  // if set to false we print out all valid/invalid barcodes.
  bool sample_bool = false; 
};

InputOptionsFastqProcess readOptionsFastqProcess(int argc, char** argv);

INPUT_OPTIONS_FASTQ_READ_STRUCTURE readOptionsFastqSlideseq(int argc, char** argv);

INPUT_OPTIONS_FASTQ_READ_STRUCTURE readOptionsFastqMetrics(int argc, char** argv);

int64_t get_num_blocks(InputOptionsFastqProcess const& options);
int64_t get_num_blocks(INPUT_OPTIONS_FASTQ_READ_STRUCTURE const& options);

#endif // __SCTOOLS_FASTQPREPROCESSING_INPUT_OPTIONS_H_
