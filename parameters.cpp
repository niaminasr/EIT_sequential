/* parameters.cpp
 * Parameters class. Input file reader.
 */
#include "parameters.hpp"
#include "stringTools.hpp"
#include <iostream>
#include <sstream>
#include <fstream>

using std::string;
using std::vector;

// ======================================================================
// Class Parameters
// ======================================================================
// PUBLIC MEMBERS
// ======================================================================

// ======================================================================
// Construct, assign

Parameters::Parameters() {}

Parameters::Parameters(string fileName) {
  _fileName = fileName;
  readInputFile();
}

Parameters::Parameters(const Parameters& that) {
  this->_fileName   = that._fileName;
  // We do not use readInputFile as the config 
  // may have been changed manually...
  this->_parameters = that._parameters;
}

Parameters& Parameters::operator=(Parameters& that) {
  if (this!=&that) {
    _fileName   = that._fileName;
    _parameters = that._parameters;
  }
  return *this;
}

// ======================================================================
// Get values from configuration

double Parameters::getReal(string key) {
  string KEY = toUpper(stringTrim(key));
  if (!hasKey(KEY))
    errNoKey(key);
  else
    return std::stod(_parameters.at(KEY));
  return 0.e0;
}

int Parameters::getInt(string key) {
  string KEY = toUpper(stringTrim(key));
  if (!hasKey(KEY))
    errNoKey(key);
  else
    return std::stoi(_parameters.at(KEY));
  return 0;
}

string Parameters::getString(string key) {
  string KEY = toUpper(stringTrim(key));
  if (!hasKey(KEY))
    errNoKey(key);
  else
    return _parameters.at(KEY);
  return string("");
}

double Parameters::getReal(string key, double defVal) {
  string KEY = toUpper(stringTrim(key));
  return hasKey(KEY) ? std::stod(_parameters.at(KEY)) : defVal;
}

int Parameters::getInt(string key, int defVal) {
  string KEY = toUpper(stringTrim(key));
  return hasKey(KEY) ? std::stoi(_parameters.at(KEY)) : defVal;
}

string Parameters::getString(string key, string defVal) {
  string KEY = toUpper(stringTrim(key));
  return hasKey(KEY) ? _parameters.at(KEY) : defVal;
}

// =====================================================================
// Modify configuration

void Parameters::set(string key,string value) {
  string KEY   = toUpper(stringTrim(key));
  string VALUE = stringTrim(value);
  _parameters[KEY] = VALUE; 
}

// =====================================================================
// Quick info

bool Parameters::hasKey(string key) {

  return _parameters.find(stringTrim(toUpper(key))) != _parameters.end();

}

// =====================================================================
// I/O

void Parameters::display() {

  std::cout << "Parameters from input file " << _fileName << std::endl;
  for (auto it: _parameters)
    std::cout << " - " << it.first << " : " << it.second << std::endl;

}

void Parameters::save(string fileName) {

  std::ofstream f(fileName.c_str());
  for (auto it: _parameters)
    f << it.first << " = " << it.second << std::endl;

}

// ======================================================================
// PRIVATE MEMBERS
// ======================================================================

void Parameters::readInputFile(bool clear) {

  if (clear) _parameters.clear();
 
  std::ifstream pfile(_fileName.c_str());
  if (!pfile.is_open())
    std::cerr << "Couldn't open parameters file: " << _fileName     << std::endl 
              << "           Will use default values if available." << std::endl;

  vector<string> lines;
  string         l,key,value,KEY;
  // Store all the lines of the file
  while(getline(pfile,l)) 
    lines.push_back(l);

  // Loop on lines
  for (string& line:lines) {
    // The line must have an = sign and must not be a comment
    if (line.find('=')!=line.npos and line.at(0)!='#') {
      std::stringstream ss(line);
      getline(ss,key,'=');
      getline(ss,value);

      key   = stringTrim(key);
      KEY   = toUpper(key);
      value = stringTrim(value);

      if (KEY.length()==0)
        std::cerr << "Warning: Empty key detected in configuration file" << std::endl;
      if (value.length()==0)
        std::cerr << "Warning: Empty value for key \""<<key<<"\"" << std::endl;
      if (_parameters.find(KEY) != _parameters.end())
        std::cerr << "Warning: Key \"" << key << "\" found several times in parameters file." << std::endl
                  << "         The value begin read last will be used." << std::endl;
    
      _parameters.insert(std::make_pair(KEY,value)); 

    }
  }

}

void Parameters::errNoKey(string key) const {
  std::cerr << "Error : no key \"" << key << "\" in configuration file and no default value provided." << std::endl;
  abort();
}

