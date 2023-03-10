#ifndef FASTQ_PREPROCESSING_WHITELIST_CORRECTOR_H_
#define FASTQ_PREPROCESSING_WHITELIST_CORRECTOR_H_

#include <string>
#include <vector>
#include <unordered_map>

#include "input_options.h"

// Manages the error correction (at most 1 Hamming distance) of raw barcodes to
// a list of expected barcodes read from a 10x Genomics whitelist file.
struct WhiteListCorrector
{
  // Maps from all correctable barcodes to indices into the 'whitelist' vector,
  // where the corresponding corrected barcode can be found. An index value of
  // -1 means the looked up key appears in the whitelist unmodified, so no need
  // to look it up in the vector. If a barcode is not present as a key, then it
  // is not correctable, i.e. >1 Hamming distance from all whitelist entries.
  //
  // "Ties" are broken in favor of whichever appeared latest in the whitelist.
  // E.g. if the whitelist file contains AT, AA, and TA in that order, then the
  // input AA would be corrected to TA, despite AA being an exact match.
  // This keeps the behavior identical to a previous Python implementation.
  // (In practice, whitelist entries are expected to be >1 Hamming distance
  //  from each other, so, famous last words, this bug should never happen.)
  std::unordered_map<std::string, int64_t> mutations;

  // all of the barcodes listed in the whitelist file, without any mutations.
  std::vector<std::string> whitelist;
};

// Builds WhiteListCorrector as described above, from the 10x Genomics whitelist
// file found at filepath white_list_file.
WhiteListCorrector readWhiteListFile(std::string const& white_list_file);

// Generates all possible single position mutations of 'barcode' and adds them
// to 'corrector', as described above.
//
// Returns false if 'barcode' has an unexpected character (i.e. not ACGTN)
bool addMutationsOfBarcodeToWhiteList(WhiteListCorrector& corrector, std::string barcode);

#endif // FASTQ_PREPROCESSING_WHITELIST_CORRECTOR_H_
