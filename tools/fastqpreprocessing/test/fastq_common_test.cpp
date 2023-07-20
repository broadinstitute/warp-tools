#include "../src/fastq_common.h"

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

// test to ensure that when an empty string is passed the returned vector is empty.
TEST(ParseReadStructureTest, EmptyStringTest) {
  std::string read_structure = "";
  std::vector<std::pair<char, int>> result = parseReadStructure(read_structure);
  EXPECT_TRUE(result.empty());
}

// test to ensure that its working as expected when passing A3 as input
TEST(ParseReadStructureTest, SingleElementTest) {
  std::string read_structure = "3A";
  std::vector<std::pair<char, int>> result = parseReadStructure(read_structure);
  ASSERT_EQ(result.size(), 1);
  EXPECT_EQ(result[0].first, 'A');
  EXPECT_EQ(result[0].second, 3);
}

// test with read structures: scATAC: 16C, V2: 16C10M, V3/multiome: 16C12M, Slideseq: 8C 18X 6C 9M 1X
TEST(ParseReadStructureTest, ReadStructureTest) {
  //expected outputs for each of the read strucutures 	
  std::vector<std::pair<char, int>> expected_result_atac = {{'C', 16}};
  std::vector<std::pair<char, int>> expected_result_v2 = {{'C', 16}, {'M', 10}};
  std::vector<std::pair<char, int>> expected_result_v3 = {{'C', 16}, {'M', 12}};
  std::vector<std::pair<char, int>> expected_result_slideseq = {{'C', 8}, {'X', 18}, {'C', 6}, {'M', 9}, {'X', 1}};

  // Test read_structure = "16C"
  std::string read_structure_atac = "16C";
  std::vector<std::pair<char, int>> result_atac = parseReadStructure(read_structure_atac);
  EXPECT_EQ(result_atac, expected_result_atac);

  // Test read_structure = "16C10M"
  std::string read_structure_v2 = "16C10M";
  std::vector<std::pair<char, int>> result_v2 = parseReadStructure(read_structure_v2);
  EXPECT_EQ(result_v2, expected_result_v2);

  // Test read_structure = "16C12M"
  std::string read_structure_v3 = "16C12M";
  std::vector<std::pair<char, int>> result_v3 = parseReadStructure(read_structure_v3);
  EXPECT_EQ(result_v3, expected_result_v3);

  // Test read_structure = "8C18X6C9M1X"
  std::string read_structure_slideseq = "8C18X6C9M1X";
  std::vector<std::pair<char, int>> result_slideseq = parseReadStructure(read_structure_slideseq);
  EXPECT_EQ(result_slideseq, expected_result_slideseq); 
}

// test reverseComplement with an test example
TEST(ReverseComplementTest, ReverseComplementSequenceTest) {
  // Test case: "ATCG" -> "CGAT"
  std::string sequence1 = "ATCG";
  std::string expected_result1 = "CGAT";
  std::string result1 = reverseComplement(sequence1);
  EXPECT_EQ(result1, expected_result1);
}

// test in case of not atac to make sure the barcode orientation is set properly 
TEST(MainCommonTest, BarcodeOrientation_FirstBPIfR3sEmptyTest) {
  // Define test inputs
  std::string white_list_file = "test/input_test_data/whitelist.txt";
  std::string barcode_orientation="FIRST_BP";
  int num_writer_threads = 1;
  std::string output_format = "BAM";
  std::vector<std::string> I1s;
  std::vector<std::string> R1s = {"test/input_test_data/R1_1.fastq", "test/input_test_data/R1_1.fastq"};
  std::vector<std::string> R2s = {"test/input_test_data/R2_1.fastq", "test/input_test_data/R2_1.fastq"};
  std::vector<std::string> R3s;  // Empty R3s vector
  std::string sample_id = "Sample1";
  std::vector<std::pair<char, int>> g_parsed_read_structure = {{'C', 16}};
  bool sample_bool = true;

  // Call the function under test
  mainCommon(white_list_file, barcode_orientation, num_writer_threads, output_format,
             I1s, R1s, R2s, R3s, sample_id, g_parsed_read_structure, sample_bool);

  // Perform necessary assertions
  EXPECT_EQ(barcode_orientation, "FIRST_BP");
}

