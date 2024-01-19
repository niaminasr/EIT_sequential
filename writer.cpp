#include "writer.hpp"
#include <fstream>
#include <unistd.h>
#include <sstream>
#include <iomanip>
#include "Eigen/Eigen/Eigen"
#include"Eigen/Eigen/Dense"

Writer::Writer(Parameters& p,Grid& g) {

  _params = &p;
  _grid   = &g;


  auto dir = p.getString("output dir");
  system(("mkdir -p " + dir).c_str());

  _baseName = dir + "/" + p.getString("output file name");

  _nFiles = 0;

}

void Writer::write(Eigen::VectorXd& u, Eigen::VectorXd& v, int i) {
  u.resize(_grid->getNbPoints(0)*_grid->getNbPoints(1) +_grid->getNbPoints(2) +_grid->Get_nbsigne());
  writeHeader(i);

  std::ofstream ofile(getFileName(i),std::ios_base::app);
  for (uint j=0;j<_grid->getNbPoints(1);j++){
    for (uint i=0;i<_grid->getNbPoints(0);i++){

       if( _grid->Levelset(i,j,v) <= 0){
         ofile << u[_grid->ind(i,j)-1] << std::endl;
       }else{
         ofile << 0 << std::endl;
        }
     }
    }

}

std::string Writer::getFileName(int i) {

  std::stringstream ss;
  // This 3 should normally be computed from t0, tf and Noutput,
  // so a link to the HeatProblem could be as useful as the link to the grid
  ss << _baseName << std::setfill('0') << std::setw(3) << i << ".vtk";
  return ss.str();

}

void Writer::writeHeader(int i) {
  std::ofstream of(getFileName(i).c_str());
  int M = _grid->getNbPoints(0);
  int N = _grid->getNbPoints(1);

  of << "# vtk DataFile Version 3.0" << std::endl;
  of << "Test heat equation" << std::endl;
  of << "ASCII" << std::endl;
  of << "DATASET STRUCTURED_POINTS" << std::endl;
  of << "DIMENSIONS " << M << " " << N << " 1" << std::endl;
  of << "ORIGIN 0. 0. 0." << std::endl;
  of << "SPACING " << _grid->getSpacing(0) << " " << _grid->getSpacing(1) << " 0." << std::endl;
  of << "POINT_DATA " << M*N << std::endl;
  of << "SCALARS u double" << std::endl;
  of << "LOOKUP_TABLE default" << std::endl;
}
