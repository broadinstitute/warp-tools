#include "fastq_common.h"
#include "input_options.h"

unsigned int g_barcode_length;
unsigned int g_umi_length;

std::string reverseComplement(std::string sequence)
{
  // get reverse complement of dna sequence 
  reverse(sequence.begin(), sequence.end());
  for (std::size_t i = 0; i < sequence.length(); ++i){
      switch (sequence[i]){
      case 'A':
        sequence[i] = 'T';
        break;    
      case 'C':
        sequence[i] = 'G';
        break;
      case 'G':
        sequence[i] = 'C';
        break;
      case 'T':
        sequence[i] = 'A';
        break;
        }
    }
  return sequence;
}

void fillSamRecord(SamRecord* samRecord, FastQFile* fastQFileI1,
                   FastQFile* fastQFileR1, FastQFile* fastQFileR2, FastQFile* fastQFileR3,
                   bool has_I1_file_list, bool has_R3_file_list,
                   std::string orientation)  
{
  // check the sequence names matching
  std::string sequence = std::string(fastQFileR1->myRawSequence.c_str());
  std::string quality_sequence = std::string(fastQFileR1->myQualityString.c_str());

  std::string barcode_seq, barcode_quality, umi_seq, umi_quality;

  // in the case of gex data -- not atac -- does not apply for slideseq (will need to change later on)
  if (!has_R3_file_list)
  {
      // extract the raw UMI
      umi_seq = sequence.substr(g_barcode_length, g_umi_length);
      // extract raw UMI quality string
      umi_quality = quality_sequence.substr(g_barcode_length, g_umi_length);
  }

  // extract the raw barcode and barcode quality 
  if (strcmp(orientation.c_str(), "FIRST_BP") == 0)
  {
      barcode_seq = sequence.substr(0, g_barcode_length);
      barcode_quality = quality_sequence.substr(0, g_barcode_length);
  }
  else if (strcmp(orientation.c_str(), "LAST_BP") == 0)
  {
      barcode_seq = sequence.substr(sequence.length() - g_barcode_length, sequence.length());
      barcode_quality = quality_sequence.substr(quality_sequence.length() - g_barcode_length, quality_sequence.length());
  }
  else if (strcmp(orientation.c_str(), "FIRST_BP_RC") == 0)
  {
      barcode_seq = reverseComplement(sequence).substr(0, g_barcode_length);
      reverse(quality_sequence.begin(), quality_sequence.end());
      barcode_quality = quality_sequence.substr(0, g_barcode_length);
  }
  else if (strcmp(orientation.c_str(), "LAST_BP_RC") == 0)
  {     
      std::string reverse_complement = reverseComplement(sequence);
      barcode = reverse_complement.substr(reverse_complement.length() - g_barcode_length, reverse_complement.length());
      
      reverse(quality_sequence.begin(), quality_sequence.end());
      barcode_quality = quality_sequence.substr(0, g_barcode_length);
  }
  else 
      crash(std::string("Incorrect barcode orientation format.\n"));
  
  fillSamRecordCommon(samRecord, fastQFileI1, fastQFileR1, fastQFileR2, fastQFileR3, 
                      has_I1_file_list, has_R3_file_list,
                      barcode_seq, barcode_quality, umi_seq, umi_quality);                
}

std::string barcodeGetter(SamRecord* sam)
{
  return std::string(sam->getString("CR").c_str());
}

void outputHandler(WriteQueue* cur_write_queue, SamRecord* samrec, int reader_thread_index)
{
  cur_write_queue->enqueueWrite(std::make_pair(samrec, reader_thread_index));
}

int main(int argc, char** argv)
{
  InputOptionsFastqProcess options = readOptionsFastqProcess(argc, argv);
  // number of output bam files, and one writer thread per bam file
  int num_writer_threads = get_num_blocks(options);
  // hardcoded this to 1000 in case of large files
  num_writer_threads =  (num_writer_threads > 1000) ? 1000 : num_writer_threads;

  g_barcode_length = options.barcode_length;
  g_umi_length = options.umi_length;

  mainCommon(options.white_list_file, options.barcode_orientation, 
             num_writer_threads, options.output_format,
             options.I1s, options.R1s, options.R2s, options.R3s, options.sample_id,
             fillSamRecord, barcodeGetter, outputHandler);

  return 0;
}
