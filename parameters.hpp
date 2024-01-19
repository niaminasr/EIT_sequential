/* parameters.hpp
 * Parameters class. Input file reader.
 */
#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#include <map>
#include <vector>
#include <string>

/* Manages options for simulations.
 *
 * Options are usually read from an input file, but can be set manually.
 * The input line format is the following:
 *
 * keyword = a list of options
 *
 * If the line ends with a backslash character (\),
 * options can be continued on the next line
 *
 */
class Parameters {

  public :

    Parameters();

    /* Constructor with a file name
     * fileName can be absolute or relative.
     */
    Parameters(std::string);

    /* Copy constructor
     * All options are copied from the configuration state,
     * not the input file (in case some entries were added/modified)
     */
    Parameters(const Parameters& cfg);//cfg=configuration.

    /* Assignment operator.
     * All options are copied from the configuration state,
     * not the input file (in case some entries were added/modified)
     */
    Parameters& operator=(Parameters& cfg);

    // Destructor. Nothing special to do.
   ~Parameters(){};

    // ===============================================================
    // Get values from config
    // ===============================================================

    // Return values associated to key. If key is missing from file,
    // program stops
    double      getReal   (std::string key);
    int         getInt    (std::string key);
    std::string getString (std::string key);

    // Returns value associated to key, if not in file, returns default value
    double      getReal   (std::string key, double      defVal);
    int         getInt    (std::string key, int         defVal);
    std::string getString (std::string key, std::string defVal);

    // ===============================================================
    // Modify configuration
    // ===============================================================

    // Replaces value in key, or creates entry if not existing
    void set (std::string key, std::string value);

    // ===============================================================
    // Quick information
    // ===============================================================

    // True if key is in the parameters map
    bool hasKey (std::string key);

    // ===============================================================
    // I/O
    // ===============================================================

    // Writes configuration in given file
    void save (std::string fileName);
    // Display configuration on standard output
    void display();

  private :

    /* Parses input file and store configuration.
     * clear: If true(default), previously defined options are deleted.
     *        When false, can be used to complete a configuration with another file.
     *        Previously existing entries are overwritten though.
     *
     *  The input line format is the following:
     *
     *  keyword = a list of options
     *
     *  If the line ends with a backslash character (\\),
     *  options can be continued on the next line
     *
     */
    void readInputFile(bool clear=true);

    // Prints an error message if file misses a key
    void errNoKey(std::string key) const;

    std::string                       _fileName;   // input file
    std::map<std::string,std::string> _parameters; // configuration info

};

#endif /* _PARAMETERS_HPP */
