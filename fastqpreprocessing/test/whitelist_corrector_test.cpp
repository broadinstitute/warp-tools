#include "../src/whitelist_corrector.h"

#include "gmock/gmock-matchers.h"
#include "gtest/gtest.h"

TEST(WhiteListCorrectorTest, BasicParsing)
{
  WhiteListCorrector corrector;
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "ACGTN"));
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "ACGGG"));
  EXPECT_FALSE(addMutationsOfBarcodeToWhiteList(corrector, " ACGTN"));
  EXPECT_FALSE(addMutationsOfBarcodeToWhiteList(corrector, "AC GTN"));
  EXPECT_FALSE(addMutationsOfBarcodeToWhiteList(corrector, "ACGTN "));
  EXPECT_FALSE(addMutationsOfBarcodeToWhiteList(corrector, "ACG\tTN"));
  EXPECT_FALSE(addMutationsOfBarcodeToWhiteList(corrector, "ACXT"));

  EXPECT_EQ(corrector.mutations.find("ACXT"), corrector.mutations.end());
  // means "exactly this key is present in the whitelist"
  EXPECT_EQ(corrector.mutations.find("ACGGG")->second, -1);
}

// (not intended to be used this way, but possible)
TEST(WhiteListCorrectorTest, DifferentLengths)
{
  WhiteListCorrector corrector;
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "ACGGG"));
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "AT"));
  // means "exactly this key is present in the whitelist"
  EXPECT_EQ(corrector.mutations.find("ACGGG")->second, -1);
  EXPECT_EQ(corrector.mutations.find("AT")->second, -1);
  EXPECT_EQ(corrector.mutations.find("GA"), corrector.mutations.end());
  EXPECT_EQ(corrector.mutations.find("GAA"), corrector.mutations.end());
}

// The logic adds whitelist entries one by one, overwriting previously added ones.
// This even includes overwriting an exact whitelist match with a correction.
// E.g. here, looking up AA ought to find the exact match AA, but instead
// "corrects" it to TA!!!
// HACK / BUG: This is bug-for-bug compatibility with a previous Python script.
// In practice this bug is expected to never happen.
TEST(WhiteListCorrectorTest, LastWhitelistEntryWinsTies)
{
  WhiteListCorrector corrector;
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "AT"));
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "AA"));
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "TA"));

  auto it = corrector.mutations.find("AA");
  EXPECT_EQ(corrector.whitelist[it->second], "TA");
}

TEST(WhiteListCorrectorTest, FindExactMatch)
{
  WhiteListCorrector corrector;
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "AT"));
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "AA"));

  // means "exactly this key is present in the whitelist"
  EXPECT_EQ(corrector.mutations.find("AA")->second, -1);
}

// I don't think we particularly care if this is supported or forbidden, but
// the logic happens to support it, so here is a test demonstrating that.
TEST(WhiteListCorrectorTest, DuplicatedWhitelistEntriesOk)
{
  WhiteListCorrector corrector;
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "CAAAT"));
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "AAAAT"));
  EXPECT_TRUE(addMutationsOfBarcodeToWhiteList(corrector, "CAAAT"));

  auto it = corrector.mutations.find("AAAAT");
  EXPECT_EQ(corrector.whitelist[it->second], "CAAAT");
}
