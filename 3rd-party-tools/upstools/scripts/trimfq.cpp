//
//  trimfq.cpp
//  upstoools
//

#include "trimfq.hpp"

void trimfq::run(string fastq, string nstart, string nend){
  // trim n bp from fastq sequence file
  int startn = stoi(nstart);
  int endn = stoi(nend);

  FILE * R1;
  string s1 = "zcat ";
  string s3 = s1 + fastq;
  R1 = popen(s3.c_str(), "r");
  string prefix = fastq.substr(0, fastq.length()-6); //.fq.gz
  
  FILE * fq_out1;
  s1 = "gzip - > ";
  s3 = s1 + prefix + "_trim.fq.gz";
  fq_out1 = popen(s3.c_str(), "w");
  
  char buffer[2000];
  while(fgets(buffer, sizeof(buffer), R1)){
    string line1(buffer);
    line1 = cxstring::chomp(line1);
    fqline rline1;
    rline1.read_part_record(R1, line1);
    if(rline1.seq.length() <= endn){
      rline1.write_record(fq_out1);
    }else{
      rline1.seq = rline1.seq.substr(startn-1, endn);
      rline1.qual = rline1.qual.substr(startn-1, endn);
      rline1.mark = "+";
      rline1.write_record(fq_out1);
    }
  }
 
  pclose(R1);
  pclose(fq_out1);

  return;
  
}
