//
//  trimfq.hpp
//  upstoools
//

#ifndef trimfq_hpp
#define trimfq_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "cxstring.hpp"

using namespace std;

class trimfq{
private:
  
public:
  static void run(string fastq, string nstart, string nend); // 
};

#endif 
