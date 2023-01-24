#include <algorithm>
#include <iostream>
#include <string>

#include "partial_file_merge.h"
#include "partial_sort.h"
#include "tagsort_input_options.h"

void warnIfNo_mitochondrial_gene_names_filename(INPUT_OPTIONS_TAGSORT const& options)
{
  if (options.mitochondrial_gene_names_filename.empty())
  {
    std::string msg =
"*** WARNING! You did not specify --mitochondrial_gene_names_filename.\n"
"Therefore, we fell back to selecting only genes beginning with 'mt-' (case\n"
"insensitive). Please write a list of all gene names you're interested in into\n"
"a file, and pass the filename with --mitochondrial_gene_names_filename.";
    std::cout << msg << std::endl;
    std::cerr << msg << std::endl;
  }
}

int main(int argc, char** argv)
{
  INPUT_OPTIONS_TAGSORT options = readOptionsTagsort(argc, argv);
  warnIfNo_mitochondrial_gene_names_filename(options);

  std::cout << "bam input " << options.bam_input << std::endl;
  std::cout << "temp folder " << options.temp_folder << std::endl;
  std::cout << "sorted output file " << options.sorted_output_file << std::endl;
  std::cout << "metric output file " << options.metric_output_file << std::endl;
  std::cout << "temp folder " << options.alignments_per_batch << std::endl;
  std::cout << "tags:" << std::endl;

  for (auto const& [tag, tag_order_num] : options.tag_order)
    std::cout << "\t" << tag << "\t" << tag_order_num << std::endl;

  // We first sort (by tag) several chunks of the data in parallel, written to
  // temporary files, then merge them to one big final sorted file.
  std::vector<std::string> partial_files = splitAndPartialSortToFiles(options);
  std::cout << "Merging " << partial_files.size() << " sorted files!"<< std::endl;
  int num_alignments = mergePartialFiles(options, partial_files);
  std::cout << "Processed " << num_alignments << " alignments." << std::endl;

  // we no longer need the partial files
  for (unsigned int i=0; i < partial_files.size(); i++)
    if (remove(partial_files[i].c_str()) != 0)
      std::cerr << "Warning: error deleting file " << partial_files[i] << std::endl;

  warnIfNo_mitochondrial_gene_names_filename(options);
  return 0;
}
