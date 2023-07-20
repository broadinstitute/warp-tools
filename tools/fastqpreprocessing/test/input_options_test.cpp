#include "../src/input_options.h"

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"
#include <filesystem>

// Tests if input parameters set are of correct type and size
TEST(ReadOptionsFastqProcessTest, BasicParsing)
{
  int argc = 19;
  char* argv[] = {"program", "--R1", "file1.fastq", "--R2", "file2.fastq", "--I1", "file3.fastq", "--R3", "file4.fastq",
                  "--sample-id", "sample1", "--output-format", "FASTQ", "--white-list", "whitelist.txt",
                  "--read-structure", "16C10M", "--bam-size", "1.5"};

  InputOptionsFastqProcess options = readOptionsFastqProcess(argc, argv);

  // Test if inputs are the correct data type
  ASSERT_TRUE((std::is_same_v<decltype(options.I1s), std::vector<std::string>>));
  ASSERT_TRUE((std::is_same_v<decltype(options.R1s), std::vector<std::string>>));
  ASSERT_TRUE((std::is_same_v<decltype(options.R2s), std::vector<std::string>>));
  ASSERT_TRUE((std::is_same_v<decltype(options.R3s), std::vector<std::string>>));
  ASSERT_TRUE((std::is_same_v<decltype(options.white_list_file), std::string>));
  ASSERT_TRUE((std::is_same_v<decltype(options.sample_id), std::string>));
  ASSERT_TRUE((std::is_same_v<decltype(options.output_format), std::string>));

  // Additional assertions to check sizes of R1, R2, R3 and I1
  ASSERT_GT(options.R1s.size(), 0);
  ASSERT_GT(options.R2s.size(), 0);
  ASSERT_EQ(options.R1s.size(), options.R2s.size());

  if (!options.I1s.empty()) {
	  ASSERT_EQ(options.I1s.size(), options.R1s.size());
   }
  
  if (!options.R3s.empty()) {
	  ASSERT_EQ(options.R3s.size(), options.R1s.size());
   } 

  ASSERT_EQ(options.R1s.size(), 1);
  ASSERT_EQ(options.R2s.size(), 1);
  ASSERT_EQ(options.I1s.size(), 1);
  ASSERT_EQ(options.R3s.size(), 1);

  // Tests input parameters
  ASSERT_EQ(options.sample_id, "sample1");
  ASSERT_EQ(options.output_format, "FASTQ");
  ASSERT_EQ(options.white_list_file, "whitelist.txt");
  ASSERT_EQ(options.read_structure, "16C10M");
  ASSERT_EQ(options.bam_size, 1.5);
}

// Tests that the file paths are correct and exist 
TEST(ReadOptionsFastqProcessTest, FilePathsExist) {
  InputOptionsFastqProcess options;

  // Set the file paths
  options.R1s = {"/warptools/fastqpreprocessing/test/input_test_data/R1_1.fastq"};
  options.R2s = {"/warptools/fastqpreprocessing/test/input_test_data/R2_1.fastq"};
  options.R3s = {"/warptools/fastqpreprocessing/test/input_test_data/R3_1.fastq"};
  options.I1s = {"/warptools/fastqpreprocessing/test/input_test_data/I1_1.fastq"};
  options.white_list_file = "/warptools/fastqpreprocessing/test/input_test_data/whitelist.txt";

  // Check if the file paths exist
  for (const auto& path : options.R1s) {
    ASSERT_TRUE(std::filesystem::is_regular_file(path));
  }
  for (const auto& path : options.R2s) {
    ASSERT_TRUE(std::filesystem::is_regular_file(path));
  }
  for (const auto& path : options.R3s) {
    ASSERT_TRUE(std::filesystem::is_regular_file(path));
  }
  for (const auto& path : options.I1s) {
    ASSERT_TRUE(std::filesystem::is_regular_file(path));
  }
  ASSERT_TRUE(std::filesystem::is_regular_file(options.white_list_file));
}

// Tests the barcode orientation and makes sure that the variable is set to FIRST_BP, LAST_BP, FIRST_BP_RC or LAST_BP_RC.
TEST(ReadOptionsFastqProcessTest, ValidBarcodeOrientation) {
  InputOptionsFastqProcess options;

  // Set the barcode_orientation to a valid value
  options.barcode_orientation = "FIRST_BP";

  // If we already set the barcode_orientation before, why are we testing it here?
  ASSERT_TRUE(options.barcode_orientation == "FIRST_BP" ||
              options.barcode_orientation == "LAST_BP" ||
              options.barcode_orientation == "FIRST_BP_RC" ||
              options.barcode_orientation == "LAST_BP_RC")
      << "Invalid barcode_orientation value: " << options.barcode_orientation;
}

