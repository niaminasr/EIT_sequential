#ifndef _WRITER_HPP_
#define _WRITER_HPP_

#include "grid.hpp"
#include "parameters.hpp"
#include <vector>
#include "Eigen/Eigen/Eigen"
#include"Eigen/Eigen/Dense"

// A class that manages outputs on regular cartesian grid
class Writer {

  public :
    Writer() {};
    Writer(Parameters& p, Grid& g);

    void write(Eigen::VectorXd& u, Eigen::VectorXd& v, int i);

  private:

    std::string getFileName(int i); // Returns the whole fileName from base and counter

    void writeHeader(int i);

    std::string _baseName; // Start of output files name
    int         _nFiles; // A counter used to index files

    Parameters* _params;
    Grid*       _grid;

};


#endif//_WRITER_HPP_


