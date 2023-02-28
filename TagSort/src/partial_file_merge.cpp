#include "partial_file_merge.h"

#include <fstream>
#include <functional>
#include <queue>

#include "metricgatherer.h"

constexpr int kLinesToReadInAChunk = 1000;

// Hands you lines one at a time from a file, reading in nice big efficient
// chunks whenever necessary.
class PartialFile
{
public:
  PartialFile(std::string const& filename) : file_(filename)
  {
    if (!file_)
      crash("ERROR failed to open the file " + filename);
    fillBuffer();
  }

  bool stillHasData()
  {
    if (sorted_.empty())
      fillBuffer();
    return !sorted_.empty();
  }

  // It is expected that you call stillHasData() immediately before every call
  // of takeNext(), to keep PartialFile's internal buffer populated.
  std::string takeNext()
  {
    std::string ret = sorted_.front();
    sorted_.pop();
    return ret;
  }

private:
  void fillBuffer()
  {
    int alignments_filled = 0;
    std::string line;
    while (std::getline(file_, line))
    {
      sorted_.push(line);
      if (++alignments_filled >= kLinesToReadInAChunk)
        break;
    }
  }

  std::ifstream file_;
  std::queue<std::string> sorted_;
};

// Helper class for merging several sorted lists into one, using
// std::priority_queue to keep track of which item to merge next.
class Merger
{
public:
  Merger() : heap_(greater_than_) {}

  // returns the appropriate next item, and an index into your PartialFile
  // array, from which you should giveInput() the next item into this Merger.
  // (If that PartialFile is now out of data, then you don't have to call
  //  giveInput() this time).
  std::pair<std::string, int> receiveNextOutput()
  {
    auto ret = heap_.top();
    heap_.pop();
    return ret;
  }

  void giveInput(std::string input, int src_file_ind)
  {
    heap_.emplace(input, src_file_ind);
  }

  bool empty() const { return heap_.empty(); }

private:
  // A custom comparator to tell our priority_queue how to deal with our std::pairs.
  // It is correct that this comparator is a > rather than <, because priority_queue
  // is a max-heap, and so wants to return the "largest" (i.e. having the comparator
  // seeing it as "before" no other item) first. So if you want 'A' to be returned
  // before 'B', your comparator needs to report 'B' as being "before" 'A'.
  std::function<bool(std::pair<std::string, int> const&,
                     std::pair<std::string, int> const&)>
      greater_than_ =
          [](std::pair<std::string, int> const& a,
             std::pair<std::string, int> const& b)
          {
            return a.first > b.first;
          };

  std::priority_queue<std::pair<std::string, int>,
                      std::vector<std::pair<std::string, int>>,
                      decltype(greater_than_) > heap_;
};

std::unique_ptr<MetricGatherer> maybeMakeMetricGatherer(INPUT_OPTIONS_TAGSORT const& options)
{
  if (!options.compute_metric)
    return nullptr;

  if (options.metric_type == MetricType::Cell)
  {
    return std::make_unique<CellMetricGatherer>(
        options.metric_output_file, options.gtf_file,
        options.mitochondrial_gene_names_filename);
  }
  else if (options.metric_type == MetricType::Gene)
    return std::make_unique<GeneMetricGatherer>(options.metric_output_file);
  else if (options.metric_type == MetricType::Umi)
    return std::make_unique<UmiMetricGatherer>(options.metric_output_file, getTagOrder(options));
  else
    crash("new MetricType enum value is not yet handled by MetricGatherer!");
  return nullptr;
}

// returns number of alignments processed
int mergePartialFiles(INPUT_OPTIONS_TAGSORT const& options,
                      std::vector<std::string> const& partial_files)
{
  int num_alignments = 0;

  std::unique_ptr<MetricGatherer> metric_gatherer = maybeMakeMetricGatherer(options);

  std::ofstream sorted_outfile;
  if (options.output_sorted_info)
  {
    sorted_outfile.open(options.sorted_output_file);
    if (!sorted_outfile)
      crash("Failed to open for writing " + options.sorted_output_file);
  }

  // We are going to merge the contents of each of these into a single sorted file.
  std::vector<PartialFile> sorted_partial_files;
  for (std::string const& fname : partial_files)
    sorted_partial_files.emplace_back(fname);

  // We use a heap to track which item to put into the final sorted file next:
  // we start it with one item from each input file.
  Merger merger;
  for (int i = 0; i < sorted_partial_files.size(); i++)
    if (sorted_partial_files[i].stillHasData())
      merger.giveInput(sorted_partial_files[i].takeNext(), i);

  while (!merger.empty())
  {
    // 'line' is the next item to go to the sorted file.
    auto [line, i] = merger.receiveNextOutput();

    if (options.output_sorted_info)
      sorted_outfile << line << "\n";
    if (metric_gatherer)
      metric_gatherer->ingestLine(line);
    num_alignments++;

    // Now that 'line' has been removed from the heap, it needs to be replaced
    // with another item from the same file (i) that it came from.
    // (If file i is out of data, then we don't need to replace).
    if (sorted_partial_files[i].stillHasData())
      merger.giveInput(sorted_partial_files[i].takeNext(), i);
  }

  if (metric_gatherer)
    metric_gatherer->outputMetricsLine(); // will have a final line of output waiting

  if (options.output_sorted_info)
  {
    std::cout << "Wrote "<< num_alignments << " alignments in sorted order to "
              << options.sorted_output_file << std::endl;
  }
  return num_alignments;
}
