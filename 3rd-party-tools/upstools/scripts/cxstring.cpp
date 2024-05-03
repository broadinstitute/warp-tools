//
//  cxstring.cpp
//  upstoools
//
//  Created by Chenxu Zhu on 3/26/22.
//

#include "cxstring.hpp"

vector<string> cxstring:: split(const string &s, const string &seperator){
    vector<string> result;
    typedef string::size_type string_size;
    string_size i = 0;

    while(i != s.size()){
        int flag = 0;
        while(i != s.size() && flag == 0){
            flag = 1;
            for(string_size x = 0; x < seperator.size(); ++x)
                if(s[i] == seperator[x]){
                    ++i;
                    flag = 0;
                    break;
                }
        }

        flag = 0;
        string_size j = i;
        while(j != s.size() && flag == 0){
            for(string_size x = 0; x < seperator.size(); ++x)
                if(s[j] == seperator[x]){
                    flag = 1;
                    break;
                }
            if(flag == 0)
                ++j;
        }
        if(i != j){
            result.push_back(s.substr(i, j-i));
            i = j;
        }
    }
    return result;
}



string cxstring:: chomp(string str){
    int i = 0;
    while(i < str.length()){
        if(str[i] == '\n' || str[i] == '\r'){
            str.erase(i, 1);
            --i;
        }
        ++i;
    }
    while(str[0] == ' '){
        str.erase(0, 1);
    }
    while(str[str.length()-1] == ' '){
        str.erase(str.length()-1, 1);
    }
    return str;
}



bool cxstring:: is_bam_header(const std::string& buffer){
    /*
     const std::regex pattern("@");
     string marker;
     marker = buffer.substr(0,1);
     return std::regex_match(marker, pattern);
     */
    if(buffer.substr(0, 1) == "@"){
        return 1;
    }
    else return 0;
}

int cxstring:: str2int(string str){
    stringstream tmp;
    int a;
    tmp << str;
    tmp >> a;
    return a;
}

string cxstring:: int2str(int i){
    stringstream tmp;
    string o;
    tmp << i;
    tmp >> o;
    return o;
}

string cxstring:: uc(string input){
  for(int i = 0; i < input.length(); i++){
      input[i] = toupper(input[i]);
  }
  return input;
}

string cxstring:: reverse_complmentary(string input) {
    reverse(input.begin(), input.end());
    for(int i = 0; i < input.length(); i++){
        if(input[i] == 'A'){
            input[i] = 'T';
        }
        else if(input[i] == 'C'){
            input[i] = 'G';
        }
        else if(input[i] == 'T'){
            input[i] = 'A';
        }
        else if(input[i] == 'G'){
            input[i] = 'C';
        }
    }
    return input;
}

int cxstring:: hamming_dist(string str1, string str2){
  // if(str1.length() != str2.length()){
  //   cout << "String not the same length:\t" << str1 << "\t" << str2 << endl;
  //   exit(256);
  // }
  int dist = 0;
  for(int i = 0; i < str1.length(); ++i){
    if(str1[i]!=str2[i])++dist;
  }
  return dist;
}



