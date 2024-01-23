#include "metricgatherer.h"
#include "mitochondrial_gene_selector.h"

#include <map>
#include <iostream>
#include <sstream>
#include <string>

constexpr int kMetricsFloatPrintPrecision = 10;

std::string to_nan(float x)
{
  std::stringstream s;
  s << std::setprecision(kMetricsFloatPrintPrecision) << x;
  return x==-1 ? "nan" : s.str();
}

MetricGatherer::MetricGatherer(std::string metric_output_file,
                               TagOrder tag_order,
                               std::string gtf_file,
                               std::string mitochondrial_gene_names_filename)
{
  std::cout<<"Constructor for metric gatherer called.\n";
  std::cout<<"Get gene ID position\n";
  std::istringstream tag_order_i(tagOrderToString(tag_order));
  std::string token;
  int geneid_position;
 
  // Tokenize tag_order_str and check positions
  for (int i = 0; std::getline(tag_order_i, token, ','); ++i) {
    if (token == "gene_id" && (i == 0 || i == 1 || i == 2)) {
        geneid_position = i; 
        std::cout << "'gene_id' is at position: " << geneid_position << std::endl;
        break;
    }
  }

  // get list of mitochondrial genes 
  if (gtf_file.empty())
    crash("MetricGatherer needs a non-empty gtf_file name!");
  // it's ok if mitochondrial_gene_names_filename is empty;
  // getInterestingMitochondrialGenes() has logic to handle that case.
  mitochondrial_genes_ = getInterestingMitochondrialGenes(
                                gtf_file, mitochondrial_gene_names_filename);

  metrics_csv_outfile_.open(metric_output_file);
  if (!metrics_csv_outfile_)
    crash("Failed to open for writing " + metric_output_file);
  metrics_csv_outfile_ << std::setprecision(kMetricsFloatPrintPrecision);
}

MetricGatherer::~MetricGatherer() {}

void MetricGatherer::clearCellAndGeneCommon()
{
  n_reads_ = 0;
  // noise_reads = 0; //# long polymers, N-sequences; NotImplemented
  fragment_histogram_.clear();
  molecule_histogram_.clear();

  molecule_barcode_fraction_bases_above_30_.clear();
  perfect_molecule_barcodes_ = 0;
  genomic_reads_fraction_bases_quality_above_30_.clear();
  genomic_read_quality_.clear();

  reads_mapped_exonic_ = 0;
  reads_mapped_exonic_as_ = 0;
  reads_mapped_intronic_ = 0;
  reads_mapped_intronic_as_ = 0;

  // alignment uniqueness information
  reads_mapped_uniquely_ = 0;
  reads_mapped_multiple_ = 0;
  duplicate_reads_ = 0;

  // alignment splicing information
  spliced_reads_ = 0;
}

void MetricGatherer::ingestLineCellAndGeneCommon(LineFields const& fields)
{
  n_reads_++;

  // the tags passed to this function define a molecule, this increments the counter,
  // identifying a new molecule only if a new tag combination is observed
  std::string hyphenated_tags = fields.tag_triple.first + "-" +
                                fields.tag_triple.second + "-" +
                                fields.tag_triple.third;
  molecule_histogram_[hyphenated_tags] += 1;

  molecule_barcode_fraction_bases_above_30_.update(fields.molecule_barcode_base_above_30);
  perfect_molecule_barcodes_ += fields.perfect_molecule_barcode;

  genomic_reads_fraction_bases_quality_above_30_.update(fields.genomic_reads_base_quality_above_30);
  genomic_read_quality_.update(fields.genomic_read_quality);

  // only do aligned-read-specific stuff if the read is mapped
  if (fields.reference != "*")
    parseAlignedReadFields(fields, hyphenated_tags);
}

void MetricGatherer::parseAlignedReadFields(LineFields const& fields, std::string hyphenated_tags)
{
  // get components that define a unique sequence fragment and increment the histogram
  std::string is_strand = (fields.is_strand == 1 ? "true" : "false");
  std::string ref_pos_str_tags = fields.reference + "\t" +
                                 std::to_string(fields.position) + "\t" +
                                 is_strand + "\t" + hyphenated_tags;
  fragment_histogram_[ref_pos_str_tags] += 1;

  // tag_order_str is a combination of BGU so find order of where gene_id is in gene_id,barcode,umi
  std::cout << "Fields tag triple\n";
  std::map<size_t, std::string> indexToField_TagOrder = {
      {0, fields.tag_triple.first}, 
      {1, fields.tag_triple.second}, 
      {2, fields.tag_triple.third}};

  std::cout << geneid_position << "\n" ;
  std::string gene_id = indexToField_TagOrder[geneid_position]; 
  std::cout << "gene name " << gene_id << "\n";
  
  // Check if not a mitochondrial gene
  if (!(mitochondrial_genes_.find(gene_id) != mitochondrial_genes_.end())) {
   if (fields.number_mappings == 1) {
      reads_mapped_uniquely_ += 1;
      if (fields.alignment_location == 1 || fields.alignment_location == 3)
        reads_mapped_exonic_ += 1;
      else if (fields.alignment_location == 2 || fields.alignment_location == 4)
        reads_mapped_exonic_as_ += 1;
      else if (fields.alignment_location == 5)
        reads_mapped_intronic_ += 1;
      else if (fields.alignment_location == 6)
        reads_mapped_intronic_as_ += 1; }
    else {
      reads_mapped_multiple_ += 1;  // without multi-mapping, this number is zero!
    }
  }
  else {
    std::cout<<"Check if mitochrondrial gene\n";
    std::cout<<"GENE " << std::string(fields.tag_triple.third) <<"\n";
    std::cout<<"mitochrondrial gene\n"; 
  }

  // in futher check if read maps outside window (when we add a  gene model)
  // and  create distances from terminate side (needs gene model) uniqueness
  duplicate_reads_ += fields.read_is_duplicate;
  spliced_reads_ += fields.read_spliced;
}

void MetricGatherer::outputMetricsLineCellAndGeneCommon()
{
  float reads_per_molecule = -1.0f;   // float("nan")
  if (molecule_histogram_.size() != 0)
    reads_per_molecule = n_reads_ / (float)molecule_histogram_.size();

  float reads_per_fragment = -1.0f; //float("nan")
  if (fragment_histogram_.size() != 0)
    reads_per_fragment = n_reads_ / (float)fragment_histogram_.size();

  float fragments_per_molecule = -1.0f; // float("nan")
  if (molecule_histogram_.size() != 0)
    fragments_per_molecule = fragment_histogram_.size() / (float)molecule_histogram_.size();

  int fragments_with_single_read_evidence = 0;
  for (auto const& [key, val] : fragment_histogram_)
    if (val == 1)
      fragments_with_single_read_evidence++;

  int molecules_with_single_read_evidence = 0;
  for (auto const& [key, val] : molecule_histogram_)
    if (val == 1)
      molecules_with_single_read_evidence++;

  metrics_csv_outfile_
      << prev_tag_ << ","
      << n_reads_ << ","
      << noise_reads << ","
      << perfect_molecule_barcodes_ << ","
      << reads_mapped_exonic_ << ","
      << reads_mapped_exonic_as_ << ","
      << reads_mapped_intronic_ << ","
      << reads_mapped_intronic_as_ << ","
      << reads_mapped_uniquely_ << ","
      << reads_mapped_multiple_ << ","
      << duplicate_reads_ << ","
      << spliced_reads_ << ","
      << kAntisenseReads << ","
      << molecule_barcode_fraction_bases_above_30_.getMean() << ","
      << to_nan(molecule_barcode_fraction_bases_above_30_.calculateVariance()) << ","
      << genomic_reads_fraction_bases_quality_above_30_.getMean() << ","
      << to_nan(genomic_reads_fraction_bases_quality_above_30_.calculateVariance()) << ","
      << to_nan(genomic_read_quality_.getMean()) << ","
      << to_nan(genomic_read_quality_.calculateVariance()) << ","
      << molecule_histogram_.size() << ","
      << fragment_histogram_.size() << ","
      << reads_per_molecule << ","
      << reads_per_fragment << ","
      << fragments_per_molecule << ","
      << fragments_with_single_read_evidence << ","
      << molecules_with_single_read_evidence;
}

int MetricGatherer::getGeneIdPosition() const{
    std::cout<<"GET GENE ID"<< geneid_position <<"\n";
    return geneid_position;
}

std::unordered_set<std::string> MetricGatherer::getMTgenes() const {
    return mitochondrial_genes_;
}

////////////////  CellMetricGatherer ////////////////////////

CellMetricGatherer::CellMetricGatherer(std::string metric_output_file,
                                      TagOrder tag_order,
                                      std::string gtf_file,
                                      std::string mitochondrial_gene_names_filename)
  : MetricGatherer(metric_output_file, tag_order, gtf_file, mitochondrial_gene_names_filename)
{
  std::cout<<"CELL METRIC GATHERER CONSTRUCTOR\n";
  std::cout<< getGeneIdPosition()<<"\n";
  
  // write metrics csv header
  std::string s;
  for (int i=0; i<25; i++)
    metrics_csv_outfile_ << "," << kCommonHeaders[i]; // TODO ok to start with ,?
  for (int i=0; i<11; i++)
    metrics_csv_outfile_ << "," << cell_specific_headers[i];
  metrics_csv_outfile_ << "\n";
}

bool MetricGatherer::cellAndGeneIsItTimeToOutput(std::string const& first_tag)
{
  // TODO?     is it ok that both cell and gene move on based only on the first tag?
  //           hmm sort of? it sounds like the user is expected to know to pass the
  //           tag types in the order that will make this work. so ideally it
  //           shouldn't be like that, but that would probably mean also remaking
  //           or even removing that "tag order" command line arg.
  return prev_tag_ != first_tag && !prev_tag_.empty();
}

void CellMetricGatherer::ingestLine(std::string const& str)
{
  LineFields fields(str);

  if (fields.tag_triple.first == "None")
    return; // ignore the None gene

  // One line of metrics file output per tag.
  if (cellAndGeneIsItTimeToOutput(fields.tag_triple.first))
  {
    outputMetricsLine();
    clear();
  }

  ingestLineCellAndGeneCommon(fields);

  // BEGIN cell-metric-specific stuff
  cell_barcode_fraction_bases_above_30_.update(fields.cell_barcode_base_above_30);
  perfect_cell_barcodes_ += fields.cell_barcode_perfect;

  // need to change this 
  // tag_order_str is a combination of BGU so find order of where gene_id is in gene_id,barcode,umi
  mitochondrial_genes_ = getMTgenes();
  geneid_position = getGeneIdPosition();
  
  std::cout << "Fields tag triple in ingestLine\n";
  std::map<size_t, std::string> indexToField_TagOrder = {
      {0, fields.tag_triple.first}, 
      {1, fields.tag_triple.second}, 
      {2, fields.tag_triple.third}};

  std::cout << geneid_position << "\n" ;
  std::string gene_id = indexToField_TagOrder[getGeneIdPosition()]; 
  std::cout << "gene name in ingestLine" << gene_id << "\n";

  if (fields.alignment_location == 7) {
    if (fields.number_mappings == 1)
      if (!(mitochondrial_genes_.find(gene_id) != mitochondrial_genes_.end()))
        reads_mapped_intergenic_ += 1;
    }
    else if(fields.alignment_location == 0) {
      reads_unmapped_ += 1;
    }

  genes_histogram_[std::string(fields.tag_triple.third)] += 1;
  // END cell-metric-specific stuff

  // TODO is it ok that we don't update prev_tag_ if reference==* ? since prev_tag_
  //      is involved in deciding when to write a metrics line, it seems like this
  //      could mess with metrics output, which includes the stuff that is updated
  //      in this function before the `if reference==* return`.
  if (fields.reference != "*")
    prev_tag_ = fields.tag_triple.first;
}

void CellMetricGatherer::outputMetricsLine()
{
  outputMetricsLineCellAndGeneCommon();

  // The number of genes that are observed by more than one read in this cell
  int genes_detected_multiple_observations = 0;
  // The number of mitochondrial genes detected by this cell
  int n_mitochondrial_genes = 0;
  // The number of molecules from mitochondrial genes detected for this cell
  int n_mitochondrial_molecules = 0;
  float pct_mitochondrial_molecules = 0.0f;

  for (auto const& [gene, count] : genes_histogram_)
  {
    if (count > 1)
      genes_detected_multiple_observations++;
    if (mitochondrial_genes_.find(gene) != mitochondrial_genes_.end())
    {
      n_mitochondrial_genes++;
      n_mitochondrial_molecules += count;
    }
  }

  if (n_mitochondrial_molecules > 0)
  {
    int tot_molecules = 0;
    for (auto const& [gene, count] : genes_histogram_)
      tot_molecules += count;

    pct_mitochondrial_molecules = 100.0f * (n_mitochondrial_molecules/(float)tot_molecules);
  }

  metrics_csv_outfile_
      << "," << perfect_cell_barcodes_
      << "," << reads_mapped_intergenic_
      << "," << reads_unmapped_
      << "," << kReadsMappedTooManyLoci
      << "," << to_nan(cell_barcode_fraction_bases_above_30_.calculateVariance())
      << "," << cell_barcode_fraction_bases_above_30_.getMean()
      << "," << genes_histogram_.size()
      << "," << genes_detected_multiple_observations
      << "," << n_mitochondrial_genes
      << "," << n_mitochondrial_molecules
      << "," << pct_mitochondrial_molecules
      << "\n";
}

void CellMetricGatherer::clear()
{
  clearCellAndGeneCommon();

  cell_barcode_fraction_bases_above_30_.clear();
  perfect_cell_barcodes_ = 0;
  reads_mapped_intergenic_ = 0;
  reads_unmapped_ = 0;
  genes_histogram_.clear();
}


////////////////  GeneMetricGatherer ////////////////////////
GeneMetricGatherer::GeneMetricGatherer(std::string metric_output_file,
                                       TagOrder tag_order,
                                       std::string gtf_file,
                                       std::string mitochondrial_gene_names_filename)
  : MetricGatherer(metric_output_file, tag_order, gtf_file, mitochondrial_gene_names_filename)
{
  
  std::cout<<"Constructor for gene metric gatherer called.\n";
  // write metrics csv header
  std::string s;
  for (int i=0; i<25; i++)
    metrics_csv_outfile_ << "," << kCommonHeaders[i]; // TODO ok to start with ,?
  for (int i=0; i<2; i++)
    metrics_csv_outfile_ << "," << gene_specific_headers[i];
  metrics_csv_outfile_ << "\n";
}

void GeneMetricGatherer::ingestLine(std::string const& str)
{
  LineFields fields(str);

  if (fields.tag_triple.first == "None")
    return; // ignore the None gene
  if (fields.tag_triple.first.find(',') != std::string::npos)
    return;

  // One line of metrics file output per tag.
  if (cellAndGeneIsItTimeToOutput(fields.tag_triple.first))
  {
    outputMetricsLine();
    clear();
  }

  ingestLineCellAndGeneCommon(fields);

  // BEGIN gene-metric-specific stuff
  cells_histogram_[std::string(fields.tag_triple.second)] += 1;
  // END gene-metric-specific stuff

  // TODO is it ok that we don't update prev_tag_ if reference==* ? since prev_tag_
  //      is involved in deciding when to write a metrics line, it seems like this
  //      could mess with metrics output, which includes the stuff that is updated
  //      in this function before the `if reference==* return`.
  if (fields.reference != "*")
    prev_tag_ = fields.tag_triple.first;
}

void GeneMetricGatherer::outputMetricsLine()
{
  outputMetricsLineCellAndGeneCommon();

  int number_cells_detected_multiple = 0;
  for (auto const& [cell, count] : cells_histogram_)
    if (count > 1)
      number_cells_detected_multiple++;

  metrics_csv_outfile_ <<  ","  << number_cells_detected_multiple
                       <<  ","  << cells_histogram_.size()
                       << "\n";
}

void GeneMetricGatherer::clear()
{
  clearCellAndGeneCommon();
  cells_histogram_.clear();
}



////////////////  UmiMetricGatherer ////////////////////////
UmiMetricGatherer::UmiMetricGatherer(std::string metric_output_file,
                                     TagOrder tag_order,
                                     std::string gtf_file,
                                     std::string mitochondrial_gene_names_filename)
  : MetricGatherer(metric_output_file, tag_order, gtf_file, mitochondrial_gene_names_filename)
{
  metrics_csv_outfile_ << tagOrderToString(tag_order) << ",count\n";
}

void UmiMetricGatherer::outputMetricsLine()
{
  metrics_csv_outfile_ << cur_histogram_triple_ << "," << cur_histogram_count_ << "\n";
}

void UmiMetricGatherer::ingestLine(std::string const& str)
{
  LineFields fields(str);

  std::string comma_separated_tags = fields.tag_triple.first + "," +
                                     fields.tag_triple.second + "," +
                                     fields.tag_triple.third;

  if (comma_separated_tags == cur_histogram_triple_)
    cur_histogram_count_++;
  else
  {
    outputMetricsLine();
    cur_histogram_triple_ = comma_separated_tags;
    cur_histogram_count_ = 1;
  }
}

void UmiMetricGatherer::clear()
{
  // (clear is not relevant, since UmiMetricGatherer::ingestLine handles its own state reset,
  //  and needs the value of comma_separated_tags to do so)
}
