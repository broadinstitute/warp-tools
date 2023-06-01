#include "whitelist_corrector.h"

#include <fstream>

// Returns false if barcode has an unexpected character (i.e. not ACGTN)
bool addMutationsOfBarcodeToWhiteList(WhiteListCorrector& corrector, std::string barcode)
{
  if (barcode.find_first_not_of("ACGTN") != std::string::npos)
    return false;

  corrector.whitelist.push_back(barcode);
  for (unsigned int i=0; i < barcode.size(); i++)
  {
    char saved = barcode[i];
    int target_whitelist_ind = corrector.whitelist.size() - 1;

    // If the mutation we're writing is already present, we just overwrite
    // what was there with the current.
    // This is done to have the same values for corrected barcodes
    // as in the python implementation.
    barcode[i] = 'A'; corrector.mutations[barcode] = target_whitelist_ind;
    barcode[i] = 'C'; corrector.mutations[barcode] = target_whitelist_ind;
    barcode[i] = 'G'; corrector.mutations[barcode] = target_whitelist_ind;
    barcode[i] = 'T'; corrector.mutations[barcode] = target_whitelist_ind;
    barcode[i] = 'N'; corrector.mutations[barcode] = target_whitelist_ind;

    barcode[i] = saved;
  }

  // -1 suggests it is already a whitelisted barcode
  // This is used, instead of the actual index, because when
  // the barcode is seen with -1 then no correction is necessary.
  // Avoids lots of vector lookups, as most barcodes are not erroneous.
  corrector.mutations[barcode] = -1;

  return true;
}

WhiteListCorrector readWhiteListFile(std::string const& white_list_file)
{
  std::ifstream file(white_list_file);
  if (!file.is_open())
    crash("Couldn't open whitelist file " + white_list_file);

  WhiteListCorrector corrector;
  for (std::string barcode; getline(file, barcode); )
    if (!addMutationsOfBarcodeToWhiteList(corrector, barcode))
      crash("Character other than ACGTN in whitelist file "+white_list_file+" line: '"+barcode+"'");

  return corrector;
}
