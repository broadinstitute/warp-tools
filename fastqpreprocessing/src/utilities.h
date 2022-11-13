#ifndef __OPTIMUS_UTILITES__
#define __OPTIMUS_UTILITES__

#include <string>
#include <vector>
#include <unordered_map>

// structure for correcting the barcodes
struct WhiteListData
{
  // Maps from barcodes to indices into the 'barcodes' vector.
  // The keys will be all barcodes in the whitelist, as well as all barcodes
  // that are Hamming distance 1 from a barcode in the whitelist.
  std::unordered_map<std::string, int64_t> mutations;
  // vector of (original) whitelist barcodes
  std::vector<std::string> barcodes;
};

// Builds barcode correction map whitelist barcodes & mutations
//
// A barcode is computed by checking if it is either in the whitelist
// (in the 10x Genomics whitelist file found at filepath white_list_file) or
// 1 mutation away from a whitelisted barcode.
WhiteListData readWhiteList(std::string const& white_list_file);

void crashWithPerror(std::string msg);
void crash(std::string msg);

#endif
