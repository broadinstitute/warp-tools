#include "partial_sort.h"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <htslib/sam.h>

#include "alignment_datatype.h"

constexpr int kThreshold = 30; // qual score threshold

extern "C" {
  bam_hdr_t* sam_hdr_read(samFile*);   //read header
  htsFile* hts_open(const char* fn, const char* mode);
}

// returns int tag, or -1 if not present
inline int get_itag_or_default(bam1_t* aln, const char* tagname, int default_value)
{
  uint8_t* p;
  int tag_value = -1;
  if ((p = bam_aux_get(aln, tagname)) == nullptr)
    tag_value = default_value;
  else
    tag_value = bam_aux2i(p);

  return  tag_value;
}

// returns string tag, or "-" if not present
inline char* get_Ztag_or_default(bam1_t* aln, const char* tagname, char* default_value)
{
  uint8_t* p;
  char* tag_value = nullptr;
  if ((p = bam_aux_get(aln, tagname)) == nullptr)
    tag_value = default_value;
  else
  {
    tag_value = bam_aux2Z(p);
    if (strcmp(tag_value, "-") == 0)
      tag_value = default_value;
  }
  return  tag_value;
}

std::unique_ptr<LineFields> parseOneAlignment(
    bam1_t* aln, INPUT_OPTIONS_TAGSORT& options, const bam_hdr_t* bam_hdr,
    TagOrder tag_order)
{
  // "consts" that the library doesn't allow to be const.
  char empty[] = "";
  char none[] = "None";
  char nochr[] = "*";

  // extract the barcodes corrected and  corrected
  char* barcode = get_Ztag_or_default(aln, options.barcode_tag.c_str(), none);
  char* barcode_raw = get_Ztag_or_default(aln, "CR", empty);

  // to be called perfect, the corrected and raw barcodes should match
  int perfect_cell_barcode = (strcmp(barcode, barcode_raw) == 0) ? 1 : 0;

  // barcode quality score
  char* barcode_qual = get_Ztag_or_default(aln, "CY", empty);

  //average barcode across the query and the fraction of barcodes above threshold
  float sum_barcode_qual = 0;
  float num_bp_above_threshold = 0;
  size_t len = strlen(barcode_qual);
  for (unsigned int k = 0; k < len; k++)
  {
    // barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
    uint8_t qual_score = (((uint8_t)barcode_qual[k]) - 33);
    sum_barcode_qual += qual_score;
    if (qual_score > kThreshold)
      num_bp_above_threshold += 1;
  }
  float avg_cell_barcode_qual = sum_barcode_qual / (float)len;
  float cell_barcode_qual_above_threshold = (float)num_bp_above_threshold / (float)len;

  // corrected molecule barcodes (UMIs)
  char* umi = get_Ztag_or_default(aln, options.umi_tag.c_str(), none);
  // raw molecule barcodes
  char* umi_raw = get_Ztag_or_default(aln, "UR", empty);

  // to be called perfect, the corrected and raw molecular barcodes should match
  int perfect_molecule_barcode = (strcmp(umi, umi_raw) == 0) ? 1 : 0;

  // qual score for molecular barcodes
  char* umi_qual = get_Ztag_or_default(aln, "UY", empty);

  float sum_umi_qual = 0;
  float num_umi_above_threshold = 0;
  len = strlen(umi_qual);
  for (unsigned int k = 0; k < len; k++)
  {
    // molecular barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
    sum_umi_qual += ((uint8_t)umi_qual[k] -33);
    if (((uint8_t)umi_qual[k] - 33) > kThreshold)
      num_umi_above_threshold += 1;
  }
  float frac_umi_qual_above_threshold = (float)num_umi_above_threshold / (float)len;

  char* gene_id = get_Ztag_or_default(aln, options.gene_tag.c_str(), none);
  char* location_tag = get_Ztag_or_default(aln, "XF", empty);

  int nh_num = get_itag_or_default(aln, "NH", -1);

  const char* chr = (aln->core.tid == -1) ? nochr : bam_hdr->target_name[aln->core.tid];

  uint32_t pos = aln->core.pos; // position.
  uint32_t isrev = bam_is_rev(aln) ? 1 : 0;   // is reverse stand
  uint32_t is_duplicate = ((aln->core.flag & BAM_FDUP) != 0) ? 1 : 0;

  // sequence quality score
  float avg_sequence_qual = 0, sum_qual = 0;
  float qual_above_threshold = 0;
  uint8_t* qual_seq = bam_get_qual(aln);  // pointer to the qual data
  len = aln->core.l_qseq; //length of qual seq.
  for (unsigned int k = 0; k < len; k++)
  {
    // the qual string are already in phred scores
    sum_qual += qual_seq[k];
    if (qual_seq[k] > kThreshold)
      qual_above_threshold += 1;
  }
  avg_sequence_qual = sum_qual / (float)len;
  qual_above_threshold = qual_above_threshold / (float)len;

  uint32_t* cigar = bam_get_cigar(aln);
  // see if it is spliced, i.e., N appears in the CIGAR string
  uint32_t spliced_read = 0;
  for (unsigned int k = 0; k < aln->core.n_cigar; k++)
  {
    uint32_t op = cigar[k] & BAM_CIGAR_MASK;
    if (op == 3 && (cigar[k] >> BAM_CIGAR_SHIFT) != 0)
    {
      spliced_read = 1;
      break;
    }
  }

  return std::make_unique<LineFields>(
      makeTriplet(barcode, umi, gene_id, tag_order), chr, location_tag, pos,
      isrev, avg_cell_barcode_qual, cell_barcode_qual_above_threshold,
      avg_sequence_qual, qual_above_threshold, nh_num, perfect_molecule_barcode,
      spliced_read, is_duplicate, perfect_cell_barcode,
      frac_umi_qual_above_threshold);
}

// Generates a random alphanumeric string (AZaz09) of a fixed length.
constexpr int kStringLen = 40;
std::string randomString()
{
  auto randchar = []() -> char
  {
    const char charset[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ rand() % max_index ];
  };
  std::string str(kStringLen, 0);
  std::generate_n(str.begin(), kStringLen, randchar);
  return str;
}


// Manages worker threads' access to reading the input file. Any worker can take
// any line from the file, but only one thread can be reading at a time. This
// class lets them take turns: they call readAlignments() whenever they want
// more data, and it blocks until they can have it.
class AlignmentReader
{
public:
  explicit AlignmentReader(INPUT_OPTIONS_TAGSORT options) : options_(options)
  {
    if ((sam_file_ptr_ = hts_open(options.bam_input.c_str(),"r")) == nullptr)
      crash(options.bam_input + ": cannot open file.");

    bam_hdr_ = sam_hdr_read(sam_file_ptr_); //read header

    allocateAlignmentBuffers();
  }

  ~AlignmentReader()
  {
    for (unsigned int i = 0; i < options_.nthreads; i++)
    {
      for (unsigned int k = 0; k < options_.alignments_per_batch; k++)
        bam_destroy1(aln_arr_[i][k]);
      free(aln_arr_[i]);
    }
    free(aln_arr_);
    sam_hdr_destroy(bam_hdr_);
    hts_close(sam_file_ptr_);
  }

  // Blocks until it's this thread's turn to read, then reads a batch of alignments.
  // Returns a pointer to the alignment ptr array, and number of alignment ptrs in the array.
  std::pair<bam1_t**, unsigned int> readAlignments(int thread_index)
  {
    const std::lock_guard<std::mutex> lock(mutex_);
    unsigned int cur_num_read = 0;
    while (cur_num_read < options_.alignments_per_batch)
    {
      // The documentation claims 0 means success, but I see code using >= and even >,
      // so it looks like actually positive means success. I guess >= is safest, then.
      if (sam_read1(sam_file_ptr_, bam_hdr_, aln_arr_[thread_index][cur_num_read]) >= 0)
        cur_num_read++;
      else
        break;
    }
    total_aligns_read_ += cur_num_read;
    batches_read_++;
    std::cout << "Finished reading batch number: " << batches_read_ << std::endl;
    return std::make_pair(aln_arr_[thread_index], cur_num_read);
  }

  void addToPartialFilenames(std::vector<std::string> names)
  {
    const std::lock_guard<std::mutex> lock(mutex_);
    for (std::string name : names)
      partial_filenames_.push_back(name);
  }

  std::vector<std::string> partial_filenames() const { return partial_filenames_; }
  bam_hdr_t* bam_hdr() const { return bam_hdr_; }
  uint64_t total_aligns_read() const { return total_aligns_read_; }

private:
  void allocateAlignmentBuffers()
  {
    assert(options_.nthreads <= kMaxTagsortThreads);
    std::string msg = "Now allocating alignment buffers. If the program crashes "
                      "here, it probably ran out of memory...";
    std::cout << msg << std::endl;
    std::cerr << msg << std::endl;

    aln_arr_ = (bam1_t***)malloc(sizeof(bam1_t**) * options_.nthreads);
    for (unsigned int i = 0; i < options_.nthreads; i++)
    {
      aln_arr_[i] = (bam1_t**)malloc(sizeof(bam1_t*) * options_.alignments_per_batch);
      for (unsigned int k = 0; k < options_.alignments_per_batch; k++)
        aln_arr_[i][k] = bam_init1(); //initialize an alignment
    }
    std::string done_msg = "Successfully allocated alignment buffers.";
    std::cout << done_msg << std::endl;
    std::cerr << done_msg << std::endl;
  }

  std::mutex mutex_;
  uint64_t total_aligns_read_ = 0;
  uint64_t batches_read_ = 0;
  INPUT_OPTIONS_TAGSORT options_;
  samFile* sam_file_ptr_ = nullptr;
  bam_hdr_t* bam_hdr_ = nullptr;
  bam1_t*** aln_arr_ = nullptr;
  std::vector<std::string> partial_filenames_;
};

void partialSortWorkerThread(int my_thread_index, AlignmentReader* alignment_reader,
                             TagOrder tag_order, INPUT_OPTIONS_TAGSORT options)
{
  std::vector<std::string> my_partial_filenames;
  bam_hdr_t* bam_hdr = alignment_reader->bam_hdr();
  while (true)
  {
    std::vector<std::unique_ptr<LineFields>> alignment_batch;

    auto [aln_ptr_array, alns_length] = alignment_reader->readAlignments(my_thread_index);
    if (alns_length == 0)
      break;

    for (unsigned int i = 0; i < alns_length; i++)
    {
      alignment_batch.push_back(parseOneAlignment(aln_ptr_array[i], options,
                                                  bam_hdr, tag_order));
    }

    // Sort our tagged alignments by their tags, and write them to a randomly
    // named text file.
    std::string partialfile_name = options.temp_folder + "/" + randomString() + ".txt";
    std::ofstream outfile(partialfile_name);

    std::sort(alignment_batch.begin(), alignment_batch.end(),
              sortAlignmentsByTagTriple);
    for (int i = 0; i < alignment_batch.size(); i++)
      alignment_batch[i]->writeTabbedToFile(outfile);

    my_partial_filenames.push_back(partialfile_name);
  }
  alignment_reader->addToPartialFilenames(my_partial_filenames);
}

std::vector<std::string> splitAndPartialSortToFiles(INPUT_OPTIONS_TAGSORT options)
{
  std::cout << "Running htslib" << std::endl;
  AlignmentReader alignment_reader(options);

  TagOrder tag_order = getTagOrder(options);

  std::vector<std::thread> worker_threads;
  for (int i = 0; i < options.nthreads; i++)
  {
    worker_threads.emplace_back(partialSortWorkerThread,
                                i, &alignment_reader, tag_order, options);
  }
  for (auto& worker_thread : worker_threads)
    worker_thread.join();

  std::cout << "Read " << alignment_reader.total_aligns_read() << " records in batches of "
            << options.alignments_per_batch << std::endl;

  return alignment_reader.partial_filenames();
}
