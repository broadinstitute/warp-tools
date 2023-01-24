#include "alignment_datatype.h"

#include <cassert>
#include <iostream>

TagOrder getTagOrder(INPUT_OPTIONS_TAGSORT options)
{
  assert(options.tag_order.size() == 3);
  // the order of the three tags are defined by the order of the supplied input arguments
  // tag.order [tag_name] -> order map
  if (options.tag_order[options.barcode_tag] == 0 &&
      options.tag_order[options.gene_tag] == 1 &&
      options.tag_order[options.umi_tag] == 2)
  {
    return TagOrder::BGU;
  }
  if (options.tag_order[options.umi_tag] == 0 &&
      options.tag_order[options.barcode_tag] == 1 &&
      options.tag_order[options.gene_tag] == 2)
  {
    return TagOrder::UBG;
  }
  if (options.tag_order[options.umi_tag] == 0 &&
      options.tag_order[options.gene_tag] == 1 &&
      options.tag_order[options.barcode_tag] == 2)
  {
    return TagOrder::UGB;
  }
  if (options.tag_order[options.gene_tag] == 0 &&
      options.tag_order[options.umi_tag] == 1 &&
      options.tag_order[options.barcode_tag] == 2)
  {
    return TagOrder::GUB;
  }
  if (options.tag_order[options.gene_tag] == 0 &&
      options.tag_order[options.barcode_tag] == 1 &&
      options.tag_order[options.umi_tag] == 2)
  {
    return TagOrder::GBU;
  }
  return TagOrder::BUG;
}

std::string tagOrderToString(TagOrder tag_order)
{
  switch (tag_order)
  {
    case TagOrder::BUG: return "barcode,umi,gene_id";
    case TagOrder::BGU: return "barcode,gene_id,umi";
    case TagOrder::UBG: return "umi,barcode,gene_id";
    case TagOrder::UGB: return "umi,gene_id,barcode";
    case TagOrder::GUB: return "gene_id,umi,barcode";
    case TagOrder::GBU: return "gene_id,barcode,umi";
    default: crash("no such TagOrder"); return "";
  }
}

TagTriple makeTriplet(std::string barcode, std::string umi, std::string gene_id,
                      TagOrder tag_order)
{
  switch (tag_order)
  {
    case TagOrder::BUG: return TagTriple(barcode, umi, gene_id);
    case TagOrder::BGU: return TagTriple(barcode, gene_id, umi);
    case TagOrder::UBG: return TagTriple(umi, barcode, gene_id);
    case TagOrder::UGB: return TagTriple(umi, gene_id, barcode);
    case TagOrder::GUB: return TagTriple(gene_id, umi, barcode);
    case TagOrder::GBU: return TagTriple(gene_id, barcode, umi);
    default: crash("no such TagOrder"); return TagTriple("","","");
  }
}

LineFields::LineFields(
    TagTriple _tag_triple, std::string _reference, std::string _alignment_location,
    int _position, int _is_strand,
    float _barcode_qual, float _cell_barcode_base_above_30,
    float _genomic_read_quality, float _genomic_reads_base_quality_above_30,
    int _number_mappings, int _perfect_molecule_barcode, int _read_spliced,
    int _read_is_duplicate, int _cell_barcode_perfect,
    float _molecule_barcode_base_above_30)
: tag_triple(_tag_triple),
  reference(_reference),
  alignment_location(_alignment_location),
  position(_position),
  is_strand(_is_strand),
  barcode_qual(_barcode_qual),
  cell_barcode_base_above_30(_cell_barcode_base_above_30),
  genomic_read_quality(_genomic_read_quality),
  genomic_reads_base_quality_above_30(_genomic_reads_base_quality_above_30),
  number_mappings(_number_mappings),
  perfect_molecule_barcode(_perfect_molecule_barcode),
  read_spliced(_read_spliced),
  read_is_duplicate(_read_is_duplicate),
  cell_barcode_perfect(_cell_barcode_perfect),
  molecule_barcode_base_above_30(_molecule_barcode_base_above_30) {}

LineFields::LineFields(std::string const& s)
{
  std::string first_tag =                           getNextField(s); // 0
  std::string second_tag =                          getNextField(s); // 1
  std::string third_tag =                           getNextField(s); // 2
  tag_triple = TagTriple(first_tag, second_tag, third_tag);
  reference =                           getNextField(s); // 3
  alignment_location =                  getNextField(s); // 4
  position =                            std::stoi(getNextField(s)); // 5
  is_strand =                           std::stoi(getNextField(s)); // 6
  barcode_qual =                        std::stof(getNextField(s)); // 7 unused
  cell_barcode_base_above_30 =          std::stof(getNextField(s)); // 8
  genomic_read_quality =                std::stof(getNextField(s)); // 9
  genomic_reads_base_quality_above_30 = std::stof(getNextField(s)); // 10
  number_mappings =                     std::stoi(getNextField(s)); // 11
  perfect_molecule_barcode =            std::stoi(getNextField(s)); // 12
  read_spliced =                        std::stoi(getNextField(s)); // 13
  read_is_duplicate =                   std::stoi(getNextField(s)); // 14
  cell_barcode_perfect =                std::stoi(getNextField(s)); // 15
  molecule_barcode_base_above_30 =      std::stof(getNextField(s)); // 16
  if (cur_tab_ != std::string::npos && cur_start_ < s.size())
    crash("Found more than the expected 17 fields in line. The bad line:\n" + s);
}

void LineFields::writeTabbedToFile(std::ofstream& outfile)
{
  outfile << tag_triple.first << "\t"
          << tag_triple.second << "\t"
          << tag_triple.third << "\t"
          << reference << "\t"
          << alignment_location << "\t"
          << position << "\t"
          << is_strand << "\t"
          << barcode_qual << "\t"
          << cell_barcode_base_above_30 << "\t"
          << genomic_read_quality << "\t"
          << genomic_reads_base_quality_above_30  << "\t"
          << number_mappings  << "\t"
          << perfect_molecule_barcode << "\t"
          << read_spliced << "\t"
          << read_is_duplicate << "\t"
          << cell_barcode_perfect << "\t"
          << molecule_barcode_base_above_30 << "\n";
}

std::string LineFields::getNextField(std::string const& s)
{
  cur_tab_ = s.find('\t', cur_start_);
  if (cur_tab_ == std::string::npos)
  {
    if (fields_gotten_ != 16)
    {
      crash("Found " + std::to_string(fields_gotten_+1) +
            " fields in line; expected 17. The bad line:\n" + s);
    }

    cur_tab_ = s.length();
    std::string ret(s.data() + cur_start_, cur_tab_ - cur_start_);
    fields_gotten_++;
    return ret;
  }
  else
  {
    std::string ret(s.data() + cur_start_, cur_tab_ - cur_start_);
    cur_start_ = cur_tab_ + 1;
    fields_gotten_++;
    return ret;
  }
}
