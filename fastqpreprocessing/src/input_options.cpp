#include "input_options.h"

#include <filesystem>
#include <getopt.h>
#include <cassert>
#include <iostream>
#include <cmath>

using std::string;

void crash(std::string msg)
{
  std::cout << msg << std::endl;
  std::cerr << msg << std::endl;
  exit(1);
}

int64_t filesize(string const& filename)
{
  FILE* f = fopen(filename.c_str(), "rb");

  int64_t size = 0;
  if (fseek(f, 0, SEEK_END) == 0)
    size = ftell(f);
  fclose(f);
  return size;
}

void printFileInfo(std::vector<string> const& fastqs,
                   string const& type)
{
  if (fastqs.size())
  {
    std::cout << "INFO " << type << " files:" << std::endl;
    for (unsigned int i= 0; i < fastqs.size(); i++)
    {
      if (std::filesystem::exists(fastqs[i].c_str()))
      {
        std::cout << "\t " << fastqs[i]  <<  " exists, file size "
                  <<  filesize(fastqs[i])  <<  std::endl;
      }
      else
      {
        std::cout << "ERROR " << fastqs[i] << " is missing!\n";
        std::cerr << "ERROR " << fastqs[i] << " is missing!\n";
        exit(1);
      }
    }
  }
}

int64_t get_num_blocks(std::vector<string> const& I1s,
                     std::vector<string> const& R1s,
                     std::vector<string> const& R2s, 
                     std::vector<string> const& R3s, double bam_size)
{
  assert(R1s.size() == R2s.size());

  if (!R3s.empty())
    assert(R1s.size() == R3s.size());

  double tot_size = 0;
  for (unsigned int i = 0; i < R1s.size(); i++)
  {
    assert(I1s.empty() || I1s.size() == R1s.size());
    if (!I1s.empty())
      tot_size += filesize(I1s[i]);

    std::cout << "file " << R1s[i] << " : " << filesize(R1s[i]) << " bytes" << std::endl;
    tot_size += filesize(R1s[i]);
    tot_size += filesize(R2s[i]);
    
    if (!R3s.empty())
      tot_size += filesize(R3s[i]);
  }

  const int GiB = 1024*1024*1024;
  return std::ceil((tot_size / GiB) / bam_size);
}

int64_t get_num_blocks(InputOptionsFastqProcess const& options)
{
  return get_num_blocks(options.I1s, options.R1s, options.R2s, options.R3s, options.bam_size);
}

int64_t get_num_blocks(INPUT_OPTIONS_FASTQ_READ_STRUCTURE const& options)
{
  return get_num_blocks(options.I1s, options.R1s, options.R2s, options.R3s, options.bam_size);
}

InputOptionsFastqProcess readOptionsFastqProcess(int argc, char** argv)
{
  InputOptionsFastqProcess options;
  int c;
  int i;
  bool verbose_flag = false;

  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"verbose",            no_argument,       0, 'v'},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"barcode-length",      required_argument, 0, 'b'},
    {"umi-length",          required_argument, 0, 'u'},
    {"bam-size",            required_argument, 0, 'B'},
    {"read-structure",      required_argument, 0, 'S'},
    {"sample-id",           required_argument, 0, 's'},
    {"I1",                  required_argument, 0, 'I'},
    {"R1",                  required_argument, 0, 'R'},
    {"R2",                  required_argument, 0, 'r'},
    {"R3",                  required_argument, 0, 'A'},
    {"barcode-orientation", required_argument, 0, 'O'},
    {"white-list",          required_argument, 0, 'w'},
    {"output-format",       required_argument, 0, 'F'},
    {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char* help_messages[] =
  {
    "verbose messages  ",
    "barcode length [required]",
    "UMI length [required]",
    "output BAM file in GB [optional: default 1 GB]",
    "read structure [required]",
    "sample id [required]",
    "I1 [optional]",
    "R1 [required -- File that contains the barcodes. This corresponds to R1 for v2/v3/multiome GEX/slideseq and R2 for scATAC.]",
    "R2 [required -- File that contains the reads. This corresponds to R2 for v2/v3/multiome GEX/slideseq. However, it corresponds to R1 in scATAC.]",
    "R3 [optional -- This file is needed for scATAC and corresponds to R2.]",
    "barcode-orientation [optional: default FIRST_BP. Other options include LAST_BP, FIRST_BP_RC or LAST_BP_RC.]",
    "whitelist (from cellranger) of barcodes [required]",
    "output-format : either FASTQ or BAM [required]",
  };


  /* getopt_long stores the option index here. */
  int option_index = 0;
  while ((c = getopt_long(argc, argv,
                          "b:u:B:s:I:R:r:w:F:v",
                          long_options,
                          &option_index)) !=- 1
        )
  {
    // process the option or arguments
    switch (c)
    {
    case 'v':
      verbose_flag = true;
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
      options.barcode_length = atoi(optarg);
      break;
    case 'u':
      options.umi_length = atoi(optarg);
      break;
    case 'B':
      options.bam_size = atof(optarg);
      break;
    case 's':
      options.sample_id = string(optarg);
      break;
    case 'S':
      options.read_structure = string(optarg);
      break;
    case 'I':
      options.I1s.push_back(string(optarg));
      break;
    case 'R':
      options.R1s.push_back(string(optarg));
      break;
    case 'r':
      options.R2s.push_back(string(optarg));
      break;
    case 'A':
      options.R3s.push_back(string(optarg));
      break;
    case 'O':
      options.barcode_orientation = string(optarg);
      break;
    case 'w':
      options.white_list_file = string(optarg);
      break;
    case 'F':
      options.output_format = string(optarg);
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
      return options;
    default:
      abort();
    }
  }

  if ((options.R1s.size() != options.R2s.size()))
  {
    crash("ERROR: Unequal number of R1 and R2 fastq files in input: R1: " +
          std::to_string(options.R1s.size()) + ", R2: " + std::to_string(options.R2s.size()));
  }

  if (options.R1s.empty())
    crash("ERROR: No R1 file provided");

  if (options.I1s.size() != options.R1s.size() && !options.I1s.empty())
    crash("ERROR: Must provide as many I1 input files as R1 input files, or else no I1 input files at all.");

  if (options.R3s.size() != options.R1s.size() && !options.R3s.empty())
    crash("ERROR: Must provide as many R3 input files as R1 input files.");

  if (options.bam_size <= 0)
    crash("ERROR: Size of a bam file (in GB) cannot be negative or 0");

  if (options.sample_id.empty())
    crash("ERROR: Must provide a sample id or name");

  if (options.output_format!="FASTQ" && options.output_format!="BAM")
    crash("ERROR: output-format must be either FASTQ or BAM");

  if (options.barcode_length <= 0)
    crash("ERROR: Barcode length must be a positive integer");

  if (options.umi_length < 0)
    crash("ERROR: UMI length must be a positive integer");

  if (verbose_flag)
  {
    if (!options.I1s.empty())
      printFileInfo(options.I1s, string("I1"));
    if (!options.R1s.empty())
      printFileInfo(options.R1s, string("R1"));
    if (!options.R2s.empty())
      printFileInfo(options.R2s, string("R2"));
    if (!options.R3s.empty())
      printFileInfo(options.R3s, string("R3"));  
  }

  return options;
}

INPUT_OPTIONS_FASTQ_READ_STRUCTURE readOptionsFastqSlideseq(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options;
  int c;
  int i;
  bool verbose_flag = false;

  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"verbose",             no_argument,       0, 'v'},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"bam-size",            required_argument, 0, 'B'},
    {"read-structure",      required_argument, 0, 'S'},
    {"sample-id",           required_argument, 0, 's'},
    {"I1",                  required_argument, 0, 'I'},
    {"R1",                  required_argument, 0, 'R'},
    {"R2",                  required_argument, 0, 'r'},
    {"R3",                  required_argument, 0, 'A'},
    {"barcode-orientation", required_argument, 0, 'O'},
    {"white-list",          required_argument, 0, 'w'},
    {"output-format",       required_argument, 0, 'F'},
    {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char* help_messages[] =
  {
    "verbose messages  ",
    "output BAM file in GB [optional: default 1 GB]",
    "read structure [required]",
    "sample id [required]",
    "I1 [optional]",
    "R1 [required -- File that contains the barcodes. This corresponds to R1 for v2/v3/multiome GEX/slideseq and R2 for scATAC.]",
    "R2 [required -- File that contains the reads. This corresponds to R2 for v2/v3/multiome GEX/slideseq. However, it corresponds to R1 in scATAC.]",
    "R3 [optional -- This file is needed for scATAC and corresponds to R2.]", 
    "barcode-orientation [optional: default FIRST_BP. Other options include LAST_BP, FIRST_BP_RC or LAST_BP_RC.]",
    "whitelist (from cellranger) of barcodes [required]",
    "output-format : either FASTQ or BAM [required]",
  };


  /* getopt_long stores the option index here. */
  int option_index = 0;
  while ((c = getopt_long(argc, argv,
                          "B:S:s:I:R:r:w:F:v",
                          long_options,
                          &option_index)) !=- 1)
  {
    // process the option or arguments
    switch (c)
    {
    case 'v':
      verbose_flag = true;
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
    case 'B':
      options.bam_size = atof(optarg);
      break;
    case 'S':
      options.read_structure = string(optarg);
      break;
    case 's':
      options.sample_id = string(optarg);
      break;
    case 'I':
      options.I1s.push_back(string(optarg));
      break;
    case 'R':
      options.R1s.push_back(string(optarg));
      break;
    case 'r':
      options.R2s.push_back(string(optarg));
      break;
    case 'A':
      options.R3s.push_back(string(optarg));
      break;
    case 'O':
      options.barcode_orientation = string(optarg);
      break;
    case 'w':
      options.white_list_file = string(optarg);
      break;
    case 'F':
      options.output_format = string(optarg);
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
      return options;
    default:
      abort();
    }
  }

  if ((options.R1s.size() != options.R2s.size()))
  {
    crash("ERROR: Unequal number of R1 and R2 fastq files in input: R1: " +
          std::to_string(options.R1s.size()) + ", R2: " + std::to_string(options.R2s.size()));
  }

  if (options.R1s.empty())
    crash("ERROR: No R1 file provided");

  if (options.I1s.size() != options.R1s.size() && !options.I1s.empty())
    crash("ERROR: Must provide as many I1 input files as R1 input files, or else no I1 input files at all.");

  if (options.R3s.size() != options.R1s.size() && !options.R3s.empty())
    crash("ERROR: Must provide as many R3 input files as R1 input files.");

  if (options.bam_size <= 0)
    crash("ERROR: Size of a bam file (in GB) cannot be negative or 0");

  if (options.sample_id.empty())
    crash("ERROR: Must provide a sample id or name");

  if (options.output_format!="FASTQ" && options.output_format!="BAM")
    crash("ERROR: output-format must be either FASTQ or BAM");

  if (options.read_structure.empty())
    crash("ERROR: Must provide read structures");

  if (verbose_flag)
  {
    if (!options.R1s.empty())
      printFileInfo(options.R1s, string("R1"));
    if (!options.R2s.empty())
      printFileInfo(options.R2s, string("R2"));
    if (!options.R3s.empty())
      printFileInfo(options.R3s, string("R3"));  
  }

  return options;
}


INPUT_OPTIONS_FASTQ_READ_STRUCTURE readOptionsFastqMetrics(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options;
  int c;
  int i;
  bool verbose_flag = false;

  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"verbose",           no_argument,       0, 'v'},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"read-structure",    required_argument, 0, 'S'},
    {"sample-id",         required_argument, 0, 's'},
    {"R1",                required_argument, 0, 'R'},
    {"white-list",        required_argument, 0, 'w'},
    {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char* help_messages[] =
  {
    "verbose messages  ",
    "read structure [required]",
    "sample id [required]",
    "R1 [required]",
    "whitelist of cell/bead barcodes [required]",
  };


  /* getopt_long stores the option index here. */
  int option_index = 0;
  while ((c = getopt_long(argc, argv,
                          "S:s:R:w:v",
                          long_options,
                          &option_index)) !=- 1
        )
  {
    // process the option or arguments
    switch (c)
    {
    case 'v':
      verbose_flag = 1;
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
    case 'S':
      options.read_structure = string(optarg);
      break;
    case 's':
      options.sample_id = string(optarg);
      break;
    case 'R':
      options.R1s.push_back(string(optarg));
      break;
    case 'w':
      options.white_list_file = string(optarg);
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
      return options;
    default:
      abort();
    }
  }

  if (options.R1s.empty())
    crash("ERROR: No R1 file provided");

  if (options.read_structure.empty())
    crash("ERROR: Must provide read structures");

  if (options.sample_id.empty())
    crash("ERROR: Must provide a sample id or name");

  if (verbose_flag && !options.R1s.empty())
    printFileInfo(options.R1s, string("R1"));

  return options;
}
