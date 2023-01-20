#ifndef TAGSORT_INPUT_OPTIONS_H_
#define TAGSORT_INPUT_OPTIONS_H_

#include <string>
#include <unordered_map>

constexpr unsigned int kMaxTagsortThreads = 30;
constexpr unsigned int kDefaultNumAlignsPerThread = 1000000;

void crash(std::string msg);

// Structure to hold input options for tagsort
enum class MetricType { Cell, Gene };
struct INPUT_OPTIONS_TAGSORT
{
  MetricType metric_type;
  bool output_sorted_info = false;
  bool compute_metric = false;
  // name of the bam file
  std::string bam_input;
  // name of the gtf file
  std::string gtf_file;
  // temp folder for disk sorting
  std::string temp_folder = "/tmp/";

  std::string metric_output_file;
  // sorted tsv output file
  std::string sorted_output_file;

  // Size (in number of alignments) of individual chunks to sort in a batch and
  // write to a partial file. Approximately 20 million alignments makes 1 GB bam file.
  unsigned int alignments_per_batch = kDefaultNumAlignsPerThread;
  unsigned int nthreads = 1;
  std::string barcode_tag;
  std::string umi_tag;
  std::string gene_tag;

  // order of the tags to sort by
  std::unordered_map<std::string, unsigned int> tag_order;

  std::string mitochondrial_gene_names_filename;
};

INPUT_OPTIONS_TAGSORT readOptionsTagsort(int argc, char** argv);

#endif // TAGSORT_INPUT_OPTIONS_H_
