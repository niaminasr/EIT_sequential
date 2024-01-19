/* stringTools.hpp
 * Utilities for std::string manipulation
 */
#ifndef _STRINGTOOLS_HPP_
#define _STRINGTOOLS_HPP_

#include <string>
#include <vector>

/* Splits s using characters from delim.
 *
 *  All trailing and preceeding spaces are removed from the original string,
 *  prior to any operation.
 *  Successive delimiters lead to empty strings being inserted.
 *  Example: stringSplit("   foo ,,, bar ,,  ",",") returns
 *  the std::vector {"foo","","","bar","",""}
 */
std::vector<std::string> stringSplit (std::string s,std::string delim=" \t");

// Removes trailing and preceeding spaces from given string.
std::string stringTrim (std::string);

// Converts eligible characters of input std::string to upper case.
std::string toUpper (std::string);

// replace
bool replace(std::string& str, const std::string& from, const std::string& to) ;
#endif /* _STRINGTOOLS_HPP_ */
