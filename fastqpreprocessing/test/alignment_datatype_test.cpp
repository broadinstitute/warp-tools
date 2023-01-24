#include "../src/alignment_datatype.h"

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(AlignmentDatatypeTest, BasicParsing)
{
  LineFields lf("a\tb\tc\td\te\t123\t1\t0.12\t4.56e10\t0\t0\t0\t0\t0\t0\t0\t0");

  EXPECT_EQ(lf.tag_triple.first, "a");
  EXPECT_EQ(lf.tag_triple.second, "b");
  EXPECT_EQ(lf.tag_triple.third, "c");
  EXPECT_EQ(lf.reference, "d");
  EXPECT_EQ(lf.alignment_location, "e");
  EXPECT_EQ(lf.position, 123);
  EXPECT_EQ(lf.is_strand, 1);
  EXPECT_NEAR(lf.barcode_qual, 0.12f, 0.001f);
  EXPECT_NEAR(lf.cell_barcode_base_above_30, 4.56e10f, 1.0f);
  EXPECT_EQ(lf.genomic_read_quality, 0);
  EXPECT_EQ(lf.genomic_reads_base_quality_above_30, 0);
  EXPECT_EQ(lf.number_mappings, 0);
  EXPECT_EQ(lf.perfect_molecule_barcode, 0);
  EXPECT_EQ(lf.read_spliced, 0);
  EXPECT_EQ(lf.read_is_duplicate, 0);
  EXPECT_EQ(lf.cell_barcode_perfect, 0);
  EXPECT_EQ(lf.molecule_barcode_base_above_30, 0);
}

TEST(AlignmentDatatypeDeathTest, FailedParseNotValidFloat)
{
  EXPECT_EXIT(LineFields lf("a\tb\tc\td\te\t1\t2\tgibberish\t1.23"),
              testing::ExitedWithCode(1),
              testing::HasSubstr("std::stof threw an exception"));
}

TEST(AlignmentDatatypeDeathTest, FailedParseCharsGivenToInt)
{
  EXPECT_EXIT(LineFields lf("a\tb\tc\td\te\tBAD123\t2"),
              testing::ExitedWithCode(1),
              testing::HasSubstr("std::stoi threw an exception"));
}

TEST(AlignmentDatatypeDeathTest, FailedParseCharsGivenAfterInt)
{
  EXPECT_EXIT(LineFields lf("a\tb\tc\td\te\t123BAD\t2"),
              testing::ExitedWithCode(1),
              testing::HasSubstr("Not an int"));
}

TEST(AlignmentDatatypeDeathTest, FailedParseFloatGivenToInt)
{
  EXPECT_EXIT(LineFields lf("a\tb\tc\td\te\t1.23\t2"),
              testing::ExitedWithCode(1),
              testing::HasSubstr("Not an int"));
}

TEST(AlignmentDatatypeDeathTest, FailedParseCharsAfterFloat)
{
  EXPECT_EXIT(LineFields lf("a\tb\tc\td\te\t1\t2\t1.23e5f"),
              testing::ExitedWithCode(1),
              testing::HasSubstr("Not a float"));
}

TEST(AlignmentDatatypeDeathTest, FailedParseTooFewFields)
{
  EXPECT_EXIT(LineFields lf("a\tb\tc\td\te\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"),
              testing::ExitedWithCode(1),
              testing::HasSubstr("Found 16 fields in line; expected 17"));
}

TEST(AlignmentDatatypeDeathTest, FailedParseTooManyFields)
{
  EXPECT_EXIT(LineFields lf("a\tb\tc\td\te\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"),
              testing::ExitedWithCode(1),
              testing::HasSubstr("Found more than the expected 17 fields in line."));
}

// TODO or is this allowed?
TEST(AlignmentDatatypeDeathTest, FailedParseEmptyField)
{
  EXPECT_EXIT(LineFields lf("a\t\tc\td\te\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"),
              testing::ExitedWithCode(1),
              testing::HasSubstr("TODO"));
}

TEST(TagTripleTest, TagOrderFromCommandLine)
{
  // TODO
//   EXPECT_EQ(getTagOrder(readOptionsTagsort(int argc, char** argv)), TagOrder::BGU);
//
//   EXPECT_EQ(getTagOrder(readOptionsTagsort(int argc, char** argv)), TagOrder::BUG);
//
//   EXPECT_EQ(getTagOrder(readOptionsTagsort(int argc, char** argv)), TagOrder::GBU);
//
//   EXPECT_EQ(getTagOrder(readOptionsTagsort(int argc, char** argv)), TagOrder::GUB);
//
//   EXPECT_EQ(getTagOrder(readOptionsTagsort(int argc, char** argv)), TagOrder::UBG);
//
//   EXPECT_EQ(getTagOrder(readOptionsTagsort(int argc, char** argv)), TagOrder::UGB);
}

TEST(TagTripleTest, MakeTagTriple)
{
  TagTriple triplet_bgu = makeTriplet("barcode", "umi", "geneID", TagOrder::BGU);
  EXPECT_EQ(triplet_bgu.first, "barcode");
  EXPECT_EQ(triplet_bgu.second, "geneID");
  EXPECT_EQ(triplet_bgu.third, "umi");

  TagTriple triplet_bug = makeTriplet("barcode", "umi", "geneID", TagOrder::BUG);
  EXPECT_EQ(triplet_bug.first, "barcode");
  EXPECT_EQ(triplet_bug.second, "umi");
  EXPECT_EQ(triplet_bug.third, "geneID");

  TagTriple triplet_gbu = makeTriplet("barcode", "umi", "geneID", TagOrder::GBU);
  EXPECT_EQ(triplet_gbu.first, "geneID");
  EXPECT_EQ(triplet_gbu.second, "barcode");
  EXPECT_EQ(triplet_gbu.third, "umi");

  TagTriple triplet_gub = makeTriplet("barcode", "umi", "geneID", TagOrder::GUB);
  EXPECT_EQ(triplet_gub.first, "geneID");
  EXPECT_EQ(triplet_gub.second, "umi");
  EXPECT_EQ(triplet_gub.third, "barcode");

  TagTriple triplet_ubg = makeTriplet("barcode", "umi", "geneID", TagOrder::UBG);
  EXPECT_EQ(triplet_ubg.first, "umi");
  EXPECT_EQ(triplet_ubg.second, "barcode");
  EXPECT_EQ(triplet_ubg.third, "geneID");

  TagTriple triplet_ugb = makeTriplet("barcode", "umi", "geneID", TagOrder::UGB);
  EXPECT_EQ(triplet_ugb.first, "umi");
  EXPECT_EQ(triplet_ugb.second, "geneID");
  EXPECT_EQ(triplet_ugb.third, "barcode");
}
