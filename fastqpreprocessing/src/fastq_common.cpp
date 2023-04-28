#include <gzstream.h>
#include <iostream>
#include <fstream>
#include <cstdint>

#include "fastq_common.h"
// number of samrecords per buffer in each reader
constexpr size_t kSamRecordBufferSize = 10000;
#include "input_options.h"
#include "whitelist_corrector.h"

#include "FastQFile.h"
#include "FastQStatus.h"
#include "BaseAsciiMap.h"
#include "SamFile.h"
#include "SamValidation.h"

#include <thread>
#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <vector>
#include <functional>
#include <stack>

// Overview of multithreading:
// * There are reader threads and writer threads. (Writers are either fastq or
//   bam, depending on how the program was run).
// * Each {reader, writer} has its own {input, output} file.
// * Each reader has an entry in g_read_arenas, and each writer has an entry in
//   g_write_queues.
// * Readers load each chunk of their processed results into SamRecord pointers
//   loaned out by their arena. They put the pointer in the correct write queue.
// * When a write queue finishes writing a SamRecord to the file, it notifies
//   the record pointer's arena that the record's memory is no longer in use.
//   The arena can then give that pointer to its reader for a new read.

PendingWrite WriteQueue::dequeueWrite()
{
  std::unique_lock<std::mutex> lock(mutex_);
  cv_.wait(lock, [&] { return !queue_.empty(); });
  auto pair = queue_.front();
  queue_.pop();
  return pair;
}
void WriteQueue::enqueueWrite(PendingWrite write)
{
  mutex_.lock();
  queue_.push(write);
  mutex_.unlock();
  cv_.notify_one();
}
void WriteQueue::enqueueShutdownSignal()
{
  mutex_.lock();
  queue_.push(std::make_pair(nullptr, kShutdown));
  mutex_.unlock();
  cv_.notify_one();
}
std::vector<std::unique_ptr<WriteQueue>> g_write_queues;

// I wrote this class to stay close to the performance characteristics of the
// original code, but I suspect the large buffers might not be necessary.
// If it doesn't slow things down noticeably, it would be cleaner to just delete
// this class, and have the WriteQueue accept unique_ptr<SamRecord> (with the
// addition of some reasonable bound on how much WriteQueue can have
// outstanding; maybe kSamRecordBufferSize items), and let them be directly
// destroyed after writing rather than be reused with this arena approach.
class SamRecordArena
{
public:
  SamRecordArena()
  {
    for (int i = 0; i < kSamRecordBufferSize; i++)
      samrecords_memory_.push_back(std::make_unique<SamRecord>());

    for (int i = samrecords_memory_.size() - 1; i >= 0; i--)
      available_samrecords_.push(samrecords_memory_[i].get());
  }

  SamRecord* acquireSamRecordMemory()
  {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [&] { return !available_samrecords_.empty(); });
    SamRecord* sam = available_samrecords_.top();
    available_samrecords_.pop();
    return sam;
  }
  void releaseSamRecordMemory(SamRecord* sam)
  {
    mutex_.lock();
    available_samrecords_.push(sam);
    mutex_.unlock();
    cv_.notify_one();
  }
private:
  std::vector<std::unique_ptr<SamRecord>> samrecords_memory_;
  std::mutex mutex_;
  std::condition_variable cv_;
  // Reusing most-recently-used memory first ought to be more cache friendly.
  std::stack<SamRecord*> available_samrecords_;
};

std::vector<std::unique_ptr<SamRecordArena>> g_read_arenas;
void releaseReaderThreadMemory(int reader_thread_index, SamRecord* samRecord)
{
  g_read_arenas[reader_thread_index]->releaseSamRecordMemory(samRecord);
}

// ---------------------------------------------------
// Write to output BAM OR FASTQ
// ----------------------------------------------------
void writeFastqRecord(ogzstream& r1_out, ogzstream& r2_out, SamRecord* sam)
{
  r1_out << "@" << sam->getReadName() << "\n" << sam->getString("CR").c_str()
         << sam->getString("UR") << "\n+\n" << sam->getString("CY") << sam->getString("UY") << "\n";
  r2_out << "@" << sam->getReadName() << "\n" << sam->getSequence() << "\n+\n"
         << sam->getQuality() << "\n";
}

void writeFastqRecordATAC(ogzstream& r1_out, ogzstream& r2_out, ogzstream& r3_out, 
                          SamRecord* sam)
{
  std::string cb_barcode = sam->getString("CB").c_str();
  std::string cr_barcode = sam->getString("CR").c_str();

  std::string write_cb_barcode = "";
  if (!cb_barcode.empty())
    write_cb_barcode = ":CB_" + cb_barcode; 
  
  //R1
  r2_out << "@" << sam->getReadName() << ":CR_" << cr_barcode << write_cb_barcode
          << "\n" << cr_barcode
          << sam->getString("UR") << "\n+\n" << sam->getString("CY") << sam->getString("UY") << "\n";
  //R2
  r1_out << "@" << sam->getReadName() << ":CR_" << cr_barcode << write_cb_barcode
          << "\n" << sam->getSequence() << "\n+\n"
          << sam->getQuality() << "\n";
  //R3
  r3_out << "@" << sam->getReadName() << ":CR_" << cr_barcode << write_cb_barcode
          << "\n" << sam->getString("RS").c_str() << "\n+\n"
          << sam->getString("RQ").c_str() <<  "\n";
  
}

void fastqWriterThread(int write_thread_index)
{
  std::string r1_output_fname = "fastq_R1_" + std::to_string(write_thread_index) + ".fastq.gz";
  ogzstream r1_out(r1_output_fname.c_str());
  if (!r1_out)
    crash("ERROR: Failed to open R1 fastq file " + r1_output_fname + " for writing");

  std::string r2_output_fname = "fastq_R2_" + std::to_string(write_thread_index) + ".fastq.gz";
  ogzstream r2_out(r2_output_fname.c_str());
  if (!r2_out)
    crash("ERROR: Failed to open R2 fastq file " + r2_output_fname + " for writing");

  while (true)
  {
    auto [sam, source_reader_index] = g_write_queues[write_thread_index]->dequeueWrite();
    if (source_reader_index == WriteQueue::kShutdown)
      break;

    writeFastqRecord(r1_out, r2_out, sam);
    g_read_arenas[source_reader_index]->releaseSamRecordMemory(sam);
  }

  // close the fastq files
  r1_out.close();
  r2_out.close();
}

//overload fastqWriterThread function for atac
void fastqWriterThreadATAC(int write_thread_index)
{
  std::string r1_output_fname = "fastq_R1_" + std::to_string(write_thread_index) + ".fastq.gz";
  ogzstream r1_out(r1_output_fname.c_str());
  if (!r1_out)
    crash("ERROR: Failed to open R1 fastq file " + r1_output_fname + " for writing");

  std::string r2_output_fname = "fastq_R2_" + std::to_string(write_thread_index) + ".fastq.gz";
  ogzstream r2_out(r2_output_fname.c_str());
  if (!r2_out)
    crash("ERROR: Failed to open R2 fastq file " + r2_output_fname + " for writing");

  std::string r3_output_fname = "fastq_R3_" + std::to_string(write_thread_index) + ".fastq.gz";
  ogzstream r3_out(r3_output_fname.c_str());
  if (!r3_out)
    crash("ERROR: Failed to open R3 fastq file " + r2_output_fname + " for writing");
  
  while (true)
  {
    auto [sam, source_reader_index] = g_write_queues[write_thread_index]->dequeueWrite();
    if (source_reader_index == WriteQueue::kShutdown)
      break;

    writeFastqRecordATAC(r1_out, r2_out, r3_out, sam);    
    g_read_arenas[source_reader_index]->releaseSamRecordMemory(sam);
  }

  // close the fastq files
  r1_out.close();
  r2_out.close();
  r3_out.close();
}

void bamWriterThread(int write_thread_index, std::string sample_id)
{
  std::string bam_out_fname = "subfile_" + std::to_string(write_thread_index) + ".bam";
  SamFile samOut;
  samOut.OpenForWrite(bam_out_fname.c_str());

  // Write the sam header.
  SamFileHeader samHeader;

  // add the HD tags for the header
  samHeader.setHDTag("VN", "1.6");
  samHeader.setHDTag("SO", "unsorted");

  // add the RG group tags
  SamHeaderRG* headerRG = new SamHeaderRG;
  headerRG->setTag("ID", "A");
  headerRG->setTag("SM", sample_id.c_str());
  samHeader.addRG(headerRG);

  // add the header to the output bam
  samOut.WriteHeader(samHeader);

  while (true)
  {
    auto [sam, source_reader_index] = g_write_queues[write_thread_index]->dequeueWrite();
    if (source_reader_index == WriteQueue::kShutdown)
      break;

    samOut.WriteRecord(samHeader, *sam);
    g_read_arenas[source_reader_index]->releaseSamRecordMemory(sam);
  }

  // close the bamfile
  samOut.Close();
}

// ---------------------------------------------------
// Parse read structure and fill sam record
// ----------------------------------------------------

// function moved from samplefastq.cpp and fastq_slideseq.cpp -- this parses the read structure
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

// get reverse complement of dna sequence -- function mainly required for atacseq data
std::string reverseComplement(std::string sequence)
{
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

// add tags to sam record 
void fillSamRecordCommon(SamRecord* samRecord, FastQFile* fastQFileI1,
                         FastQFile* fastQFileR1, FastQFile* fastQFileR2, FastQFile* fastQFileR3,
                         bool has_I1_file_list, bool has_R3_file_list,
                         std::string const& barcode_seq, std::string const& barcode_quality,
                         std::string const& umi_seq, std::string const& umi_quality)
{
  // reset the samrecord
  samRecord->resetRecord();
  // add read group and the sam flag
  samRecord->addTag("RG", 'Z', "A");
  samRecord->setFlag(4);
  // add identifier, sequence and quality score of the alignments
  samRecord->setReadName(fastQFileR2->mySequenceIdentifier.c_str());
  samRecord->setSequence(fastQFileR2->myRawSequence.c_str());
  samRecord->setQuality(fastQFileR2->myQualityString.c_str());
  // add barcode and quality
  samRecord->addTag("CR", 'Z', barcode_seq.c_str());
  samRecord->addTag("CY", 'Z', barcode_quality.c_str());
  // add UMI
  samRecord->addTag("UR", 'Z', umi_seq.c_str());
  samRecord->addTag("UY", 'Z', umi_quality.c_str());
  // add raw sequence and quality sequence for the index
  if (has_I1_file_list)
  {
    samRecord->addTag("SR", 'Z', fastQFileI1->myRawSequence.c_str());
    samRecord->addTag("SY", 'Z', fastQFileI1->myQualityString.c_str());
  }
  // add raw sequence and quality sequence for the R3 atac fastq file 
  if (has_R3_file_list)
  { 
    samRecord->addTag("RS", 'Z', fastQFileR3->myRawSequence.c_str());
    samRecord->addTag("RQ", 'Z', fastQFileR3->myQualityString.c_str());
  }
}

// fill sam record -- this function was modified and moved from fastqprocess.cpp, samplefastq.cpp and fastq_slideseq.cpp
void fillSamRecord(SamRecord* samRecord, FastQFile* fastQFileI1,
                   FastQFile* fastQFileR1, FastQFile* fastQFileR2, FastQFile* fastQFileR3,
                   bool has_I1_file_list, bool has_R3_file_list, std::string orientation, 
                   std::vector<std::pair<char, int>> g_parsed_read_structure)  
{
  // check the sequence names matching
  std::string sequence = std::string(fastQFileR1->myRawSequence.c_str());
  std::string quality_sequence = std::string(fastQFileR1->myQualityString.c_str());
  std::string barcode_seq, barcode_quality, umi_seq, umi_quality;
  int g_barcode_length;

  // extract the raw barcode and barcode quality  
  // when orientation is set to FIRST_BP use the g_parse_read_structure
  // other cases are for other atac barcode variations which depends on other factors and does not need 
  // (1) atac data (when has_R3_file_list is set to True) -- read_structure for atac will be 16C
  // (2) with slideseq/gex data (when has_R3_file_list is set to False)
  if (strcmp(orientation.c_str(), "FIRST_BP") == 0)
  {
      int cur_ind = 0;
      for (auto [tag, length] : g_parsed_read_structure)
      {
        switch (tag)
        {
          case 'C':
            barcode_seq += sequence.substr(cur_ind, length);
            barcode_quality += quality_sequence.substr(cur_ind, length);
            break;
          case 'M':
            umi_seq += sequence.substr(cur_ind, length);
            umi_quality += quality_sequence.substr(cur_ind, length);
            break;
          default:
            break;
        }
        cur_ind += length;
      }
  }
  else if (has_R3_file_list)
  {
      //with atacseq data read strucuture will look like this "16C" -- where 16 is the barcode length
      g_barcode_length = std::get<1>(g_parsed_read_structure[0]);
      
      if (strcmp(orientation.c_str(), "LAST_BP") == 0)
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
          barcode_seq = reverse_complement.substr(reverse_complement.length() - g_barcode_length, reverse_complement.length());
          
          reverse(quality_sequence.begin(), quality_sequence.end());
          barcode_quality = quality_sequence.substr(0, g_barcode_length);
      }
      else 
          crash(std::string("Incorrect barcode orientation format.\n"));
  }

  fillSamRecordCommon(samRecord, fastQFileI1, fastQFileR1, fastQFileR2, fastQFileR3, 
                      has_I1_file_list, has_R3_file_list,
                      barcode_seq, barcode_quality, umi_seq, umi_quality);                
}

// get barcode from sam record -- this function was modified and moved from fastqprocess.cpp, samplefastq.cpp and fastq_slideseq.cpp
std::string barcodeGetter(SamRecord* sam)
{
  return std::string(sam->getString("CR").c_str());
}

// ---------------------------------------------------
// Correct whitelist
// ---------------------------------------------------

// Computes the whitelist-corrected barcode and adds it to sam_record.
// Returns the index of the bamfile bucket / writer thread where sam_record
// should be sent.
int32_t correctBarcodeToWhitelist(
    const std::string& barcode, SamRecord* sam_record, const WhiteListCorrector* corrector,
    int* n_barcode_corrected, int* n_barcode_correct, int* n_barcode_errors, int num_writer_threads)
{
  std::string correct_barcode;
  // bucket barcode is used to pick the target bam file
  // This is done because in the case of incorrectible barcodes
  // we need a mechanism to uniformly distribute the alignments
  // so that no bam is oversized to putting all such barcode less
  // sequences into one particular. Incorrectible barcodes are simply
  // added without the CB tag
  std::string bucket_barcode;
  if (auto it = corrector->mutations.find(barcode) ; it != corrector->mutations.end())
  {
    int64_t mutation_index = it->second;
    if (mutation_index == -1) // -1 means raw barcode is correct
    {
      correct_barcode = barcode;
      *n_barcode_correct += 1;
    }
    else
    {
      // it is a 1-mutation of some whitelist barcode so get the
      // barcode by indexing into the vector of whitelist barcodes
      correct_barcode = corrector->whitelist[mutation_index];
      *n_barcode_corrected += 1;
    }
    // is used for computing the file index
    bucket_barcode = correct_barcode;

    // corrected barcode should be added to the samrecord
    sam_record->addTag("CB", 'Z', correct_barcode.c_str());
  }
  else     // not possible to correct the raw barcode -- aseel: is this raw?
  {
    *n_barcode_errors += 1;
    bucket_barcode = barcode;
  }
  // destination bam file index computed based on the bucket_barcode
  return std::hash<std::string> {}(bucket_barcode) % num_writer_threads;
}

// ---------------------------------------------------
// Read FASTQ one read at a time
// ---------------------------------------------------

// Returns true if successfully read a sequence.
bool readOneItem(FastQFile& fastQFileI1, bool has_I1_file_list,
                   FastQFile& fastQFileR1, FastQFile& fastQFileR2,
                   FastQFile& fastQFileR3, bool has_R3_file_list)
{
  return (!has_I1_file_list ||
      (
        has_I1_file_list &&
        fastQFileI1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS
      )
  )
  && fastQFileR1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS
  && fastQFileR2.readFastQSequence() == FastQStatus::FASTQ_SUCCESS
  && (!has_R3_file_list || 
        (
          has_R3_file_list &&
          fastQFileR3.readFastQSequence() == FastQStatus::FASTQ_SUCCESS
        )
  );
}

void fastQFileReaderThread(
    int reader_thread_index, std::string filenameI1, String filenameR1,
    String filenameR2, std::string filenameR3, const WhiteListCorrector* corrector, std::string barcode_orientation,
    std::vector<std::pair<char, int>> g_parsed_read_structure,
    std::function<void(WriteQueue*, SamRecord*, int)> output_handler)
{
  /// setting the shortest sequence allowed to be read
  FastQFile fastQFileI1(4, 4);
  FastQFile fastQFileR1(4, 4);
  FastQFile fastQFileR2(4, 4);
  FastQFile fastQFileR3(4, 4);

  bool has_I1_file_list = true;
  if (!filenameI1.empty())
  {
    if (fastQFileI1.openFile(String(filenameI1.c_str()), BaseAsciiMap::UNKNOWN) !=
        FastQStatus::FASTQ_SUCCESS)
    {
      crash(std::string("Failed to open file: ") + filenameI1);
    }
  }
  else
    has_I1_file_list = false;

  //This is for the 3rd atacseq file. 
  bool has_R3_file_list = true;
  if (!filenameR3.empty())
  {
    if (fastQFileR3.openFile(String(filenameR3.c_str()), BaseAsciiMap::UNKNOWN) !=
        FastQStatus::FASTQ_SUCCESS)
    {
      crash(std::string("Failed to open file: ") + filenameR3);
    }
  }
  else
    has_R3_file_list = false;

  if (fastQFileR1.openFile(filenameR1, BaseAsciiMap::UNKNOWN) !=
      FastQStatus::FASTQ_SUCCESS)
  {
    crash(std::string("Failed to open file: ") + filenameR1.c_str());
  }
  if (fastQFileR2.openFile(filenameR2, BaseAsciiMap::UNKNOWN) !=
      FastQStatus::FASTQ_SUCCESS)
  {
    crash(std::string("Failed to open file: ") + filenameR2.c_str());
  }

  // Keep reading the file until there are no more fastq sequences to process.
  int total_reads = 0;
  int n_barcode_errors = 0;
  int n_barcode_corrected = 0;
  int n_barcode_correct = 0;
  printf("Opening the thread in %d\n", reader_thread_index);

  while (fastQFileR1.keepReadingFile())
  {
    if (readOneItem(fastQFileI1, has_I1_file_list, fastQFileR1, fastQFileR2, fastQFileR3, has_R3_file_list))
    {
      total_reads++;

      SamRecord* samrec = g_read_arenas[reader_thread_index]->acquireSamRecordMemory();

      // prepare the samrecord with the sequence, barcode, UMI, and their quality sequences
      fillSamRecord(samrec, &fastQFileI1, &fastQFileR1, &fastQFileR2, &fastQFileR3, has_I1_file_list,
                        has_R3_file_list, barcode_orientation, g_parsed_read_structure); 

      // get barcode 
      std::string barcode = barcodeGetter(samrec);
                                    
      // bucket barcode is used to pick the target bam file
      // This is done because in the case of incorrigible barcodes
      // we need a mechanism to uniformly distribute the alignments
      // so that no bam is oversized to putting all such barcode less
      // sequences into one particular. Incorregible barcodes are simply
      // added withouth the CB tag
      int32_t bam_bucket = correctBarcodeToWhitelist(
          barcode, samrec, corrector, &n_barcode_corrected, 
          &n_barcode_correct, &n_barcode_errors, g_write_queues.size());

      output_handler(g_write_queues[bam_bucket].get(), samrec, reader_thread_index);

      if (total_reads % 10000000 == 0)
      {
        printf("%d\n", total_reads);
        std::string a = std::string(fastQFileR1.myRawSequence.c_str());
        printf("%s\n", fastQFileR1.mySequenceIdLine.c_str());
        printf("%s\n", fastQFileR2.mySequenceIdLine.c_str());
        printf("%s\n", fastQFileR3.mySequenceIdLine.c_str());
      }
    }
  }

  // Finished processing all of the sequences in the file.
  // Close the input files.
  if (has_I1_file_list)
    fastQFileI1.closeFile();
  if (has_R3_file_list)
    fastQFileR3.closeFile();  
  
  fastQFileR1.closeFile();
  fastQFileR2.closeFile();

  printf("Total barcodes:%d\n correct:%d\ncorrected:%d\nuncorrectible"
         ":%d\nuncorrected:%lf\n",
         total_reads, n_barcode_correct, n_barcode_corrected, n_barcode_errors,
         n_barcode_errors/static_cast<double>(total_reads) * 100);
}

// ---------------------------------------------------
// Main 
// ---------------------------------------------------

void mainCommon(
    std::string white_list_file, std::string barcode_orientation,
    int num_writer_threads, std::string output_format,
    std::vector<std::string> I1s, std::vector<std::string> R1s, 
    std::vector<std::string> R2s, std::vector<std::string> R3s,
    std::string sample_id,  std::vector<std::pair<char, int>> g_parsed_read_structure,
    std::function<void(WriteQueue*, SamRecord*, int)> output_handler)
{
  std::cout << "reading whitelist file " << white_list_file << "...";
  // stores barcode correction map and vector of correct barcodes
  WhiteListCorrector corrector = readWhiteListFile(white_list_file);
  std::cout << "done" << std::endl;

  for (int i = 0; i < R1s.size(); i++)
    g_read_arenas.push_back(std::make_unique<SamRecordArena>());
  for (int i = 0; i < num_writer_threads; i++)
    g_write_queues.push_back(std::make_unique<WriteQueue>());

  // execute the bam file writers threads
  std::vector<std::thread> writers;
  if (output_format == "BAM")
    for (int i = 0; i < num_writer_threads; i++)
      writers.emplace_back(bamWriterThread, i, sample_id);
  else if (output_format == "FASTQ")
    for (int i = 0; i < num_writer_threads; i++)
      if (R3s.empty())
          writers.emplace_back(fastqWriterThread, i);
      else
          writers.emplace_back(fastqWriterThreadATAC, i);
  else
    crash("ERROR: Output-format must be either FASTQ or BAM");

  // execute the fastq readers threads
  std::vector<std::thread> readers;

  for (unsigned int i = 0; i < R1s.size(); i++)
  {
    assert(I1s.empty() || I1s.size() == R1s.size());
    // if there is no I1/R3 file then send an empty file name
    readers.emplace_back(fastQFileReaderThread, i, I1s.empty() ? "" : I1s[i], R1s[i].c_str(),
                         R2s[i].c_str(), R3s.empty() ? "" : R3s[i].c_str(), 
                         &corrector, barcode_orientation,
                         g_parsed_read_structure,
                         //sam_record_filler, barcode_getter, 
                         output_handler);
  }

  for (auto& reader : readers)
    reader.join();

  // Now that there's nothing left to read, we can safely append a shutdown
  // signal to all the write queues.
  for (auto& write_queue : g_write_queues)
    write_queue->enqueueShutdownSignal();

  for (auto& writer : writers)
    writer.join();
}
