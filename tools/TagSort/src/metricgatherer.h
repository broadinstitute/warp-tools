#ifndef __METRIC_GATHERER__
#define __METRIC_GATHERER__

#include <unordered_map>
#include <string>
#include <iostream>
#include <vector>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <unordered_set>

#include "alignment_datatype.h"
#include "tagsort_input_options.h"

class OnlineGaussianSufficientStatistic
{
public:
  // incorporates new_value into the online estimate of mean and variance
  void update(double new_value)
  {
    _count += 1.0;
    _sum += new_value;
    sum_EX2 += (new_value*new_value);
  }

  // return the mean value
  double getMean()
  {
    _mean = _sum/_count;
    return _mean;
  }

  // calculate and return the variance
  double calculateVariance()
  {
    if (_count < 2)
      return -1.0;
    return sum_EX2 / (_count - 1) - (_sum/_count) * (_sum / (_count - 1));
  }

  void clear()
  {
    _mean_squared_error = 0.0;
    _mean = 0.0;
    _count = 0;
    _sum = 0;
    sum_EX2 = 0.0;
  }

private:
  double _mean_squared_error = 0.0;
  double sum_EX2 = 0.0;
  double _mean = 0.0;
  double _sum = 0.0;
  double _count = 0.0;
};

// Base class for recording metrics from a stream of (sorted-by-tag-triple)
// alignment lines.
class MetricGatherer
{
public:
  MetricGatherer(std::string metric_output_file, TagOrder tag_order,
                 std::string gtf_file,
                 std::string mitochondrial_gene_names_filename);
  virtual ~MetricGatherer();

  virtual void ingestLine(std::string const& str) = 0;
  virtual void outputMetricsLine() = 0;

protected:
  // Each line of metric output is built from all alignments with a given tag.
  // After that output is written, the gatherer resets its recording state for
  // the next group. This function does that.
  virtual void clear() = 0;
  bool cellAndGeneIsItTimeToOutput(std::string const& first_tag);
  void ingestLineCellAndGeneCommon(LineFields const& fields);
  void outputMetricsLineCellAndGeneCommon();
  void clearCellAndGeneCommon();
  bool isMitochondrial(LineFields const& fields) const;

  const std::string kCommonHeaders[27] =
  {
    "n_reads",
    "tso_reads",
    "noise_reads",
    "perfect_molecule_barcodes",
    "reads_mapped_exonic",
    "reads_mapped_exonic_as",
    "reads_mapped_intronic",
    "reads_mapped_intronic_as",
    "reads_mapped_uniquely",
    "reads_mapped_multiple",
    "duplicate_reads",
    "spliced_reads",
    "antisense_reads",
    "molecule_barcode_fraction_bases_above_30_mean",
    "molecule_barcode_fraction_bases_above_30_variance",
    "genomic_reads_fraction_bases_quality_above_30_mean",
    "genomic_reads_fraction_bases_quality_above_30_variance",
    "genomic_read_quality_mean",
    "genomic_read_quality_variance",
    "n_molecules",
    "n_fragments",
    "reads_per_molecule",
    "reads_per_fragment",
    "fragments_per_molecule",
    "fragments_with_single_read_evidence",
    "molecules_with_single_read_evidence",
    "reads_mapped_mitochondrial"
  };

  void parseAlignedReadFields(LineFields const& fields, std::string hyphenated_tags);

  std::unordered_set<std::string> mito_genes_;
  const TagOrder tag_order_;

  std::unordered_map<std::string, int> molecule_histogram_;
  std::ofstream metrics_csv_outfile_;

  std::string prev_tag_;

private:
  // count information
  int n_reads_ = 0;
  int tso_reads_ = 0;
  const int noise_reads = 0; //# long polymers, N-sequences; NotImplemented

  std::unordered_map<std::string, int> fragment_histogram_;

  // molecule information
  OnlineGaussianSufficientStatistic molecule_barcode_fraction_bases_above_30_;

  int perfect_molecule_barcodes_ = 0;

  OnlineGaussianSufficientStatistic genomic_reads_fraction_bases_quality_above_30_;

  OnlineGaussianSufficientStatistic genomic_read_quality_;

  // (Note that all of these reads_mapped fields count only unique reads; any
  //  read that has duplicates does not contribute to them *at all*)
  // alignment location information
  int reads_mapped_exonic_ = 0;
  int reads_mapped_exonic_as_ = 0;
  int reads_mapped_intronic_ = 0;
  int reads_mapped_intronic_as_ = 0;
  //int reads_mapped_utr_ = 0;
  int reads_mapped_mitochondrial_ = 0;

  // in future we can implement this when we have a gene model
  // self.reads_mapped_outside_window = 0  # reads should be within 1000 bases of UTR
  // self._read_distance_from_termination_site = OnlineGaussianSufficientStatistic()

  // alignment uniqueness information
  int reads_mapped_uniquely_ = 0;
  int reads_mapped_multiple_ = 0;
  int duplicate_reads_ = 0;

  // alignment splicing information
  int spliced_reads_ = 0;
  const int kAntisenseReads = 0; // TODO is never changed from 0
  // int plus_strand_reads_ = 0;  // strand balance (currently unused)
  // test for mt -- will need to remove
  int n_mitochondrial_reads_ = 0;
};

class CellMetricGatherer: public MetricGatherer
{
public:
  CellMetricGatherer(std::string metric_output_file,
                     TagOrder tag_order,
                     std::string gtf_file,
                     std::string mitochondrial_gene_names_filename);
  void ingestLine(std::string const& str) override;
  void outputMetricsLine() override;

protected:
  void clear() override;

private:
  int perfect_cell_barcodes_ = 0; // The number of reads whose cell barcodes contain no errors (tag ``CB`` == ``CR``)
  int reads_mapped_intergenic_ = 0; // The number of reads mapped to an intergenic region for this cell

  int reads_unmapped_ = 0;

  //  The number of reads that were mapped to too many loci across the genome and as a
  //  consequence, are reported unmapped by the aligner
  const int kReadsMappedTooManyLoci = 0; // TODO is never changed from 0

  OnlineGaussianSufficientStatistic cell_barcode_fraction_bases_above_30_;
  std::unordered_map<std::string, int> genes_histogram_;

  std::string cell_specific_headers[11] =
  {
    "perfect_cell_barcodes",
    "reads_mapped_intergenic",
    "reads_unmapped",
    "reads_mapped_too_many_loci",
    // The variance of the fraction of Illumina base calls for the cell barcode
    // sequence that are greater than 30, across molecules
    "cell_barcode_fraction_bases_above_30_variance",
    // The average fraction of Illumina base calls for the cell barcode sequence
    // that are greater than 30, across molecules
    "cell_barcode_fraction_bases_above_30_mean",
    "n_genes",
    "genes_detected_multiple_observations",
    "n_mitochondrial_genes",
    "n_mitochondrial_molecules",
    "pct_mitochondrial_molecules"
  };
};


class GeneMetricGatherer: public MetricGatherer
{
public:
  GeneMetricGatherer(std::string metric_output_file,
                     TagOrder tag_order,
                     std::string gtf_file,
                     std::string mitochondrial_gene_names_filename);

  void ingestLine(std::string const& str) override;

  void outputMetricsLine() override;

protected:
  void clear() override;

private:
  std::unordered_map<std::string, int> cells_histogram_;
  std::string gene_specific_headers[2] =
  {
    "number_cells_detected_multiple",
    "number_cells_expressing"
  };
};

class UmiMetricGatherer: public MetricGatherer
{
public:
  UmiMetricGatherer(std::string metric_output_file,
                    TagOrder tag_order,
                    std::string gtf_file,
                    std::string mitochondrial_gene_names_filename);
  void ingestLine(std::string const& str) override;
  void outputMetricsLine() override;

protected:
  void clear() override;

private:
  std::string cur_histogram_triple_{};
  int cur_histogram_count_ = 0;
};

#endif
