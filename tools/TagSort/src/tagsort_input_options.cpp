#include "tagsort_input_options.h"

#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <string>

using std::string;

void crash(string msg)
{
  std::cout << msg << std::endl;
  std::cerr << msg << std::endl;
  exit(1);
}

INPUT_OPTIONS_TAGSORT readOptionsTagsort(int argc, char** argv)
{
  INPUT_OPTIONS_TAGSORT options;
  int c;
  int i;

  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"compute-metric",             no_argument,       0, 'm'},
    {"output-sorted-info",         no_argument,       0, 'n'},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"bam-input",                  required_argument, 0, 'b'},
    {"gtf-file",                   required_argument, 0, 'a'},
    {"temp-folder",                required_argument, 0, 't'},
    {"sorted-output",              required_argument, 0, 'o'},
    {"metric-output",              required_argument, 0, 'M'},
    {"alignments-per-thread",      required_argument, 0, 'p'},
    {"nthreads",                   required_argument, 0, 'T'},
    {"barcode-tag",                required_argument, 0, 'C'},
    {"umi-tag",                    required_argument, 0, 'U'},
    {"gene-tag",                   required_argument, 0, 'G'},
    {"metric-type",                required_argument, 0, 'K'},
    {"mitochondrial-gene-names-filename", required_argument, 0, 'g'},
    {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char* help_messages[] =
  {
    "compute metric, metrics are computed if this option is provided [optional]",
    "sorted output file is produced if this option is provided [optional]",
    "input bam file [required]",
    "gtf file (unzipped) required then metric type is cell [required with metric cell]",
    "temp folder for disk sorting [options: default /tmp]",
    "sorted output file [optional]",
    "metric file, the metrics are output in this file  [optional]",
    "number of alignments per thread [optional: default 1000000], if this number is increased then more RAM is required but reduces the number of file splits",
    "number of threads [optional: default 1]",
    "barcode-tag the call barcode tag [required]",
    "umi-tag the umi tag [required]: the tsv file output is sorted according the tags in the options barcode-tag, umi-tag or gene-tag",
    "gene-tag the gene tag [required]",
    "metric type, one of \"cell\", \"gene\", or \"umi\" [required]",
    "file listing gene names, one per line, that the program should care about. [required, may omit if you want mouse or human]"
  };

  string metric_type_str;

  /* getopt_long stores the option index here. */
  int option_index = 0;
  int curr_size = 0;
  while ((c = getopt_long(argc, argv,
                          "b:a:t:no:mM:p:T:C:U:G:K:",
                          long_options,
                          &option_index)) !=- 1)
  {
    // process the option or arguments
    switch (c)
    {
    case 'm':
      options.compute_metric = true;
      break;
    case 'n':
      options.output_sorted_info = true;
      break;
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      printf("option %s", long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;
    case 'b':
      options.bam_input = string(optarg);
      break;
    case 'a':
      options.gtf_file = string(optarg);
      break;
    case 't':
      options.temp_folder = string(optarg);
      break;
    case 'o':
      options.sorted_output_file = string(optarg);
      break;
    case 'M':
      options.metric_output_file = string(optarg);
      break;
    case 'p':
      options.alignments_per_batch = atoi(optarg);
      break;
    case 'T':
      options.nthreads = atoi(optarg);
      break;
    case 'C':
      options.barcode_tag = string(optarg);
      curr_size = options.tag_order.size();
      options.tag_order[string(optarg)] = curr_size;
      break;
    case 'U':
      options.umi_tag = string(optarg);
      curr_size = options.tag_order.size();
      options.tag_order[string(optarg)] = curr_size;
      break;
    case 'G':
      options.gene_tag = string(optarg);
      curr_size = options.tag_order.size();
      options.tag_order[string(optarg)] = curr_size;
      break;
    case 'K':
      metric_type_str = string(optarg);
      break;
    case 'g':
      options.mitochondrial_gene_names_filename = string(optarg);
      break;
    case '?':
    case 'h':
      i = 0;
      printf("Usage: %s [options] \n", argv[0]);
      while (long_options[i].name != 0)
      {
        printf("\t--%-20s  %-25s  %-35s\n", long_options[i].name,
               long_options[i].has_arg == no_argument?
               "no argument" : "required_argument",
               help_messages[i]);
        i = i + 1;
      }
      /* getopt_long already printed an error message. */
      exit(0);
      break;
    default:
      abort();
    }
  }

  // Check the options
  // either metric computation or the sorted tsv file must be produced
  if (!options.output_sorted_info && !options.compute_metric)
    crash("ERROR: The choice of either the sorted alignment info or metric computation must be specified");

  if (options.compute_metric && options.metric_output_file.empty())
    crash("ERROR: Must specify --metric-output when specifying --compute-metric");

  if (options.output_sorted_info && options.sorted_output_file.empty())
    crash("ERROR: Must specify --sorted-output when specifying --output-sorted-info");

  if (metric_type_str == "cell")
    options.metric_type = MetricType::Cell;
  else if (metric_type_str == "gene")
    options.metric_type = MetricType::Gene;
  else if (metric_type_str == "umi")
    options.metric_type = MetricType::Umi;
  else
    crash("ERROR: --metric-type must be \"cell\", \"gene\", or \"umi\"");

  // if metric type is cell then the gtf file must be provided
  if (options.metric_type == MetricType::Cell && options.gtf_file.empty())
    crash("ERROR: The gtf file name must be provided with metric_type \"cell\"");

  // the gtf file should not be gzipped
  if (options.gtf_file.length() > 3)
  {
    int len = options.gtf_file.length();
    string const& fname = options.gtf_file;
    if (fname[len-3] == '.' && tolower(fname[len-2]) == 'g' && tolower(fname[len-1]) == 'z')
      crash("ERROR: The gtf file must not be gzipped");
  }

  // bam input file must be there
  if (options.bam_input.empty())
    crash("ERROR: Must specify a input file name");

  // check for input file
  if (!std::filesystem::exists(options.bam_input.c_str()))
    crash("ERROR: bam_input " + options.bam_input + " is missing!");

  // check for the temp folder
  if (!std::filesystem::exists(options.temp_folder.c_str()))
    crash("ERROR: temp folder " + options.temp_folder + " is missing!");

  // check for three distinct tags, barcode, umi and gene_id tags
  if (options.tag_order.size() != 3)
    crash("ERROR:  Must have three distinct tags");
  bool seen_tag_index[3] = { false, false, false };
  for (auto [tag, index] : options.tag_order)
  {
    if (index < 0 || index > 2)
      crash("Invalid tag index " + std::to_string(index) + "; must be 0 1 or 2");
    else
      seen_tag_index[index] = true;
  }
  if (!(seen_tag_index[0] && seen_tag_index[1] && seen_tag_index[2]))
    crash("Need tag indices 0 1 and 2");

  // The size of a set of aligments for in-memory sorting must be positive
  if (options.alignments_per_batch < 1000)
    crash("ERROR: The number of alignments per thread must be at least 1000");

  // The number of threads must be between 1 and kMaxTagsortThreads
  if (options.nthreads > kMaxTagsortThreads || options.nthreads < 1)
    crash("ERROR: The number of threads must be between 1 and " + std::to_string(kMaxTagsortThreads));

  return options;
}
