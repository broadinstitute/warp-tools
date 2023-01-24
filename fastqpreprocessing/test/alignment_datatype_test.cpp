#include "../src/alignment_datatype.h"

#include "gtest/gtest.h"

TEST(AlignmentDatatypeTest, BasicParsing)
{
  LineFields lf("TODO\tGET\tA\tGOOD\tEXAMPLE");

  // TODO
//   EXPECT_EQ(lf.tag_triple,
//   EXPECT_EQ(lf.reference,
//   EXPECT_EQ(lf.alignment_location,
//   EXPECT_EQ(lf.position,
//   EXPECT_EQ(lf.is_strand,
//   EXPECT_EQ(lf.barcode_qual,
//   EXPECT_EQ(lf.cell_barcode_base_above_30,
//   EXPECT_EQ(lf.genomic_read_quality,
//   EXPECT_EQ(lf.genomic_reads_base_quality_above_30,
//   EXPECT_EQ(lf.number_mappings,
//   EXPECT_EQ(lf.perfect_molecule_barcode,
//   EXPECT_EQ(lf.read_spliced,
//   EXPECT_EQ(lf.read_is_duplicate,
//   EXPECT_EQ(lf.cell_barcode_perfect,
//   EXPECT_EQ(lf.molecule_barcode_base_above_30,

}

TEST(AlignmentDatatypeTest, FailedParseNotValidFloat)
{
  bool thrown = false;
  try { LineFields lf("a\tb\tc\td\te\t1\t2\tgibberish\t1.23"); }
  catch (std::exception const& e) { thrown = true; }
  if (!thrown)
    FAIL() <<"Didn't throw an exception when it should have";
}

TEST(AlignmentDatatypeTest, FailedParseFloatGivenToInt)
{
  bool thrown = false;
  try { LineFields lf("a\tb\tc\td\te\t1.23\t2"); }
  catch (std::exception const& e) { thrown = true; }
  if (!thrown)
    FAIL() << "Didn't throw an exception when it should have";
}

TEST(AlignmentDatatypeDeathTest, FailedParseTooFewFields)
{
  EXPECT_EXIT(LineFields lf("a\tb\tc\td\te\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"),
              testing::ExitedWithCode(1),
              "Found 16 fields in line; expected 17");
}

TEST(AlignmentDatatypeDeathTest, FailedParseTooManyFields)
{
  EXPECT_EXIT(LineFields lf("a\tb\tc\td\te\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"),
              testing::ExitedWithCode(1),
              "Found 18 fields in line; expected 17");
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
