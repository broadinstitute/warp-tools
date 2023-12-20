//
//  sepType_DPT.cpp
//  2nt pre-indexde of DPT

#include "sepType_DPT.hpp"

void sepType_DPT::run(string prefix, string BC01){

  FILE * R1;
  FILE * R2;
  FILE * R3;
  int barcode = stoi(BC01);
  string s1 = "zcat ";
  string s3 = s1 + prefix + "_R1.fq.gz";
  R1 = popen(s3.c_str(), "r");

  s3 = s1 + prefix + "_R2.fq.gz";
  R2 = popen(s3.c_str(), "r");

  s3 = s1 + prefix + "_R3.fq.gz";
  R3 = popen(s3.c_str(), "r");

  FILE * Rout1;
  FILE * Rout2;
  FILE * Rout3;

  s1 = "gzip - > ";
  s3 = s1 + prefix + "_R1_prefix.fq.gz";
  Rout1 = popen(s3.c_str(), "w");
  s3 = s1 + prefix + "_R2_prefix.fq.gz";
  Rout2 = popen(s3.c_str(), "w");
  s3 = s1 + prefix + "_R3_prefix.fq.gz";
  Rout3 = popen(s3.c_str(), "w");
  
  char buffer[2000];
  while(fgets(buffer, sizeof(buffer), R1)){
    string line1(buffer);
    line1 = cxstring::chomp(line1);
    fqline rline1;
    fqline rline2;
    fqline rline3;
    rline1.read_part_record(R1, line1);
    rline2.read_full_record(R2);
    rline3.read_full_record(R3);
    vector<string> sp = cxstring::split(rline1.readname, " ");
    vector<string> sp1 = cxstring::split(rline2.readname, " ");
    vector<string> sp2 = cxstring::split(rline3.readname, " ");

    string bc = rline2.seq.substr(0, barcode);
    rline1.readname = sp[0] + ":" + bc;
    rline2.readname = sp1[0] + ":" + bc;
    rline2.seq = rline2.seq.substr(barcode, rline2.seq.length() - barcode);
    rline2.qual = rline2.qual.substr(barcode, rline2.qual.length() - barcode);
    rline3.readname = sp2[0] + ":" + bc;

    rline1.write_record(Rout1);
    rline2.write_record(Rout2);
    rline3.write_record(Rout3);
  }
  
  //close read files
  pclose(R1);
  pclose(R2);
  pclose(R3);
  
  //close output files
  pclose(Rout1);
  pclose(Rout2);
  pclose(Rout3);
  
  return;  
}
