/* stringTools.cpp
 * Utilities for string manipulation
 */
#include "stringTools.hpp"
#include <unistd.h>
#include <stdio.h>
#include <algorithm>

using std::vector;
using std::string;

vector<string> stringSplit(string s,string delims) {

  string s2 = stringTrim(s);
  vector<string> res;
  size_t posBegin,posEnd;
  // Build a string of non space delimiters 
  // (that can be occuring several times:
  // " data1 ,,, data4, ...
  string delims2(delims);
  delims2.erase(remove(delims2.begin(),delims2.end(),' ' ),delims2.end());
  delims2.erase(remove(delims2.begin(),delims2.end(),'\t'),delims2.end());

  // String starting with delimiters...
  if (delims.find(s2.front())!=delims.npos)
    res.push_back("");

  posBegin = s2.find_first_not_of(delims.c_str());
  //for (uint_wp i=0;i<posBegin;i++)
  //  res.push_back("");

  while(posBegin!=s2.npos) {
    posEnd = s2.find_first_of(delims.c_str(),posBegin);
    if (posEnd>posBegin)
      res.push_back(stringTrim(s2.substr(posBegin,posEnd-posBegin)));
    posBegin = s2.find_first_not_of(delims.c_str(),posEnd);
    // But there was maybe several successive non space delimiters
    if (posBegin!=s2.npos)
      for (size_t posT=posEnd+1;posT<posBegin;posT++)
        if (delims2.find(s2[posT])!=delims2.npos)
          res.push_back("");
  }

  // String ending with a delimiter...
  if (delims.find(s2.back())!=delims.npos)
    res.push_back("");

  return res;

}

string stringTrim(string s) {

  size_t posBegin = s.find_first_not_of(" \t");
  size_t posEnd   = s.find_last_not_of (" \t")+1;
  return (posEnd>posBegin) ?
    s.substr(posBegin,posEnd-posBegin) :
    string("") ;

}

string toUpper(string s) {
  string s2(s);
  for (auto &c: s2)
    c = toupper(c);
  return s2;
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}
