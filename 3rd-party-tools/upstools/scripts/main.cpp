//
//  main.cpp
//  upstoools
//

#include <iostream>
#include <string>
#include "cxstring.hpp"
#include "sepType_DPT.hpp"
#include "trimfq.hpp"

using namespace std;

void help(){
  // version 2022.03.27
  string ver = "2023.03.03";
  cout << "=== UPSTools ===" << endl;
  cout << "Version: " << ver << endl;
  cout << endl;
  cout << "Usage:" << endl;
  cout << "upstools sepType_DPT [prefix] [length]\n\tExtract 1st barcode (preindex) from the head of fastq sequence, append at readname." << endl;
  cout << "upstools trimfq [input_fq] [start_bp] [end_bp]\n\tTrim a fastq (end with .fq.gz) to desired length." << endl;
  cout << endl;
  
  return;
}

int main(int argc, const char * argv[]) {
  if(argc<3){
    help();
    return 0;
  }
  string mod(argv[1]);

  if(mod == "sepType_DPT"){
    if(argc<4){
      help();
      return 0;
    }
    string input(argv[2]);
    string lec(argv[3]);
    sepType_DPT::run(input, lec);
  }

  if(mod == "trimfq"){
    if(argc<5){
      help();
      return 0;
    }
    trimfq::run(argv[2], argv[3], argv[4]);
  }
  
  return 0;
}

