#include "../src/alignment_datatype.h"
#include "../src/partial_file_merge.h"

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(PartialFileMergeTest, BasicMerge)
{
  INPUT_OPTIONS_TAGSORT options;
  options.compute_metric = false;
  options.output_sorted_info = true;
  options.sorted_output_file = "/tmp/TagSort_gtest_partial_file_merge_test_sorted_output_file";

  options.barcode_tag = "barcode";
  options.umi_tag = "umi";
  options.gene_tag = "gene";
  options.tag_order["barcode"] = 0;
  options.tag_order["umi"] = 1;
  options.tag_order["gene"] = 2;

  options.alignments_per_batch = 2;
  options.nthreads = 2;

  std::vector<std::string> my_partial_filenames;
  {
    std::vector<std::unique_ptr<LineFields>> alignment_batch;
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("A", "A", "A");
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("B", "B", "B");

    std::sort(alignment_batch.begin(), alignment_batch.end(),
              sortAlignmentsByTagTriple);

    std::string partialfile_name = "/tmp/TagSort_gtest_partial_file_merge_test_temp_file_0";
    std::ofstream outfile(partialfile_name);
    my_partial_filenames.push_back(partialfile_name);
    for (int i = 0; i < alignment_batch.size(); i++)
      alignment_batch[i]->writeTabbedToFile(outfile);
  }
  {
    std::vector<std::unique_ptr<LineFields>> alignment_batch;
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("C", "C", "C");
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("B", "B", "B");

    std::sort(alignment_batch.begin(), alignment_batch.end(),
              sortAlignmentsByTagTriple);

    std::string partialfile_name = "/tmp/TagSort_gtest_partial_file_merge_test_temp_file_1";
    std::ofstream outfile(partialfile_name);
    my_partial_filenames.push_back(partialfile_name);
    for (int i = 0; i < alignment_batch.size(); i++)
      alignment_batch[i]->writeTabbedToFile(outfile);
  }

  int num_alignments_merged = mergePartialFiles(options, my_partial_filenames);
  EXPECT_EQ(num_alignments_merged, 4);

  std::ifstream sorted_file(options.sorted_output_file);
  std::string line;
  std::getline(sorted_file, line);
  EXPECT_EQ(line[0], 'A');
  std::getline(sorted_file, line);
  EXPECT_EQ(line[0], 'B');
  std::getline(sorted_file, line);
  EXPECT_EQ(line[0], 'B');
  std::getline(sorted_file, line);
  EXPECT_EQ(line[0], 'C');


  unlink("/tmp/TagSort_gtest_partial_file_merge_test_sorted_output_file");
  unlink("/tmp/TagSort_gtest_partial_file_merge_test_temp_file_0");
  unlink("/tmp/TagSort_gtest_partial_file_merge_test_temp_file_1");
}

TEST(PartialFileMergeTest, EmptyTagsOk)
{
  INPUT_OPTIONS_TAGSORT options;
  options.compute_metric = false;
  options.output_sorted_info = true;
  options.sorted_output_file = "/tmp/TagSort_gtest_partial_file_merge_test_sorted_output_file";

  options.barcode_tag = "barcode";
  options.umi_tag = "umi";
  options.gene_tag = "gene";
  options.tag_order["barcode"] = 0;
  options.tag_order["umi"] = 1;
  options.tag_order["gene"] = 2;

  options.alignments_per_batch = 2;
  options.nthreads = 2;

  std::vector<std::string> my_partial_filenames;
  {
    std::vector<std::unique_ptr<LineFields>> alignment_batch;
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("A", "A", "A");
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("B", "B", "B");
    alignment_batch.emplace_back(std::make_unique<LineFields>()); // empty

    std::sort(alignment_batch.begin(), alignment_batch.end(),
              sortAlignmentsByTagTriple);

    std::string partialfile_name = "/tmp/TagSort_gtest_partial_file_merge_test_temp_file_0";
    std::ofstream outfile(partialfile_name);
    my_partial_filenames.push_back(partialfile_name);
    for (int i = 0; i < alignment_batch.size(); i++)
      alignment_batch[i]->writeTabbedToFile(outfile);
  }
  {
    std::vector<std::unique_ptr<LineFields>> alignment_batch;
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("C", "C", "C");
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("B", "B", "B");

    std::sort(alignment_batch.begin(), alignment_batch.end(),
              sortAlignmentsByTagTriple);

    std::string partialfile_name = "/tmp/TagSort_gtest_partial_file_merge_test_temp_file_1";
    std::ofstream outfile(partialfile_name);
    my_partial_filenames.push_back(partialfile_name);
    for (int i = 0; i < alignment_batch.size(); i++)
      alignment_batch[i]->writeTabbedToFile(outfile);
  }

  int num_alignments_merged = mergePartialFiles(options, my_partial_filenames);
  EXPECT_EQ(num_alignments_merged, 5);

  std::ifstream sorted_file(options.sorted_output_file);
  std::string line;
  std::getline(sorted_file, line);
  EXPECT_EQ(line[0], '\t'); // empty tag, so first char is the tab delim
  std::getline(sorted_file, line);
  EXPECT_EQ(line[0], 'A');
  std::getline(sorted_file, line);
  EXPECT_EQ(line[0], 'B');
  std::getline(sorted_file, line);
  EXPECT_EQ(line[0], 'B');
  std::getline(sorted_file, line);
  EXPECT_EQ(line[0], 'C');

  unlink("/tmp/TagSort_gtest_partial_file_merge_test_sorted_output_file");
  unlink("/tmp/TagSort_gtest_partial_file_merge_test_temp_file_0");
  unlink("/tmp/TagSort_gtest_partial_file_merge_test_temp_file_1");
}

// Exercises comparing tags beyond the first field.
// Tests that the sorting gives AAA -> AAB -> ABB -> BAA -> BBB
TEST(PartialFileMergeTest, TagTiebreaking)
{
  INPUT_OPTIONS_TAGSORT options;
  options.compute_metric = false;
  options.output_sorted_info = true;
  options.sorted_output_file = "/tmp/TagSort_gtest_partial_file_merge_test_sorted_output_file";

  options.barcode_tag = "barcode";
  options.umi_tag = "umi";
  options.gene_tag = "gene";
  options.tag_order["barcode"] = 0;
  options.tag_order["umi"] = 1;
  options.tag_order["gene"] = 2;

  options.alignments_per_batch = 2;
  options.nthreads = 2;

  std::vector<std::string> my_partial_filenames;
  {
    std::vector<std::unique_ptr<LineFields>> alignment_batch;
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("A", "A", "A");
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("A", "A", "B");

    std::sort(alignment_batch.begin(), alignment_batch.end(),
              sortAlignmentsByTagTriple);

    std::string partialfile_name = "/tmp/TagSort_gtest_partial_file_merge_test_temp_file_0";
    std::ofstream outfile(partialfile_name);
    my_partial_filenames.push_back(partialfile_name);
    for (int i = 0; i < alignment_batch.size(); i++)
      alignment_batch[i]->writeTabbedToFile(outfile);
  }
  {
    std::vector<std::unique_ptr<LineFields>> alignment_batch;
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("A", "B", "B");
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("B", "A", "A");
    alignment_batch.emplace_back(std::make_unique<LineFields>());
    alignment_batch.back()->tag_triple = TagTriple("B", "B", "B");

    std::sort(alignment_batch.begin(), alignment_batch.end(),
              sortAlignmentsByTagTriple);

    std::string partialfile_name = "/tmp/TagSort_gtest_partial_file_merge_test_temp_file_1";
    std::ofstream outfile(partialfile_name);
    my_partial_filenames.push_back(partialfile_name);
    for (int i = 0; i < alignment_batch.size(); i++)
      alignment_batch[i]->writeTabbedToFile(outfile);
  }

  int num_alignments_merged = mergePartialFiles(options, my_partial_filenames);
  EXPECT_EQ(num_alignments_merged, 5);

  std::ifstream sorted_file(options.sorted_output_file);
  std::string line;
  std::getline(sorted_file, line);
  EXPECT_THAT(line, testing::StartsWith("A\tA\tA"));
  std::getline(sorted_file, line);
  EXPECT_THAT(line, testing::StartsWith("A\tA\tB"));
  std::getline(sorted_file, line);
  EXPECT_THAT(line, testing::StartsWith("A\tB\tB"));
  std::getline(sorted_file, line);
  EXPECT_THAT(line, testing::StartsWith("B\tA\tA"));
  std::getline(sorted_file, line);
  EXPECT_THAT(line, testing::StartsWith("B\tB\tB"));


  unlink("/tmp/TagSort_gtest_partial_file_merge_test_sorted_output_file");
  unlink("/tmp/TagSort_gtest_partial_file_merge_test_temp_file_0");
  unlink("/tmp/TagSort_gtest_partial_file_merge_test_temp_file_1");
}
