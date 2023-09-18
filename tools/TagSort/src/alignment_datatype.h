#ifndef TAGSORT_ALIGNMENT_DATATYPE_H_
#define TAGSORT_ALIGNMENT_DATATYPE_H_

#include <fstream>
#include <memory>
#include <string>

#include "tagsort_input_options.h"

struct TagTriple
{
public:
  TagTriple() {}
  TagTriple(std::string one, std::string two, std::string three)
    : first(one), second(two), third(three) {}
  std::string first;
  std::string second;
  std::string third;
};

// TODO better name
struct LineFields
{
public:
  // Parses a LineFields from an ASCII string. 's' should be a single line of
  // tab separated values, one for each field in this struct.
  explicit LineFields(std::string const& s);

  LineFields(
      TagTriple _tag_triple, std::string _reference, int _alignment_location,
      int _position, int _is_strand,
      float _barcode_qual, float _cell_barcode_base_above_30,
      float _genomic_read_quality, float _genomic_reads_base_quality_above_30,
      int _number_mappings, int _perfect_molecule_barcode, int _read_spliced,
      int _read_is_duplicate, int _cell_barcode_perfect,
      float _molecule_barcode_base_above_30);

  LineFields() {}

  void writeTabbedToFile(std::ofstream& outfile);

  TagTriple tag_triple; // (0,1,2) barcode umi and gene_id, not necessarily in that order
  std::string reference; // 3
  int alignment_location; // (4) aka biotype (TODO which is more accurate? or is this wrong?)
  int position; // 5
  int is_strand; // (6) 1 for yes, 0 for no
  float barcode_qual; // 7
  float cell_barcode_base_above_30; // 8
  float genomic_read_quality; // 9
  float genomic_reads_base_quality_above_30; // 10
  int number_mappings; // 11
  int perfect_molecule_barcode; // (12) 1 for yes, 0 for no
  // cigar N field (3) indicates a read is spliced if the value is non-zero
  int read_spliced; // (13) 1 for yes, 0 for no
  int read_is_duplicate; // (14) 1 for yes, 0 for no
  int cell_barcode_perfect; // (15) 1 for yes, 0 for no
  float molecule_barcode_base_above_30; // (16) fraction of umi qual score > 30
};

bool sortAlignmentsByTagTriple(std::unique_ptr<LineFields> const& a,
                               std::unique_ptr<LineFields> const& b);

// Parses tab-separated fields from a line (std::string s).
class LineFieldsParser
{
public:
  explicit LineFieldsParser(std::string const& s) : s_(s) {}
  int getNextFieldInt();
  float getNextFieldFloat();
  std::string getNextField();
  bool hasMore() const;
  void crashLF(std::string msg);

private:
  std::string const& s_;
  size_t cur_start_ = 0;
  size_t cur_tab_ = 0;
  int fields_gotten_ = 0;
};

enum class TagOrder { BUG, BGU, UBG, UGB, GUB, GBU };
TagOrder getTagOrder(INPUT_OPTIONS_TAGSORT options);
std::string tagOrderToString(TagOrder tag_order);

TagTriple makeTriplet(std::string barcode, std::string umi, std::string gene_id,
                      TagOrder tag_order);

#endif // TAGSORT_ALIGNMENT_DATATYPE_H_
