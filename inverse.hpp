#ifndef _INVERSE_HPP_
#define _INVERSE_HPP_

#include "parameters.hpp"
#include "grid.hpp"
#include "assembler.hpp"
#include "writer.hpp"
#include "heatProblem.hpp"
#include "Eigen/Eigen/IterativeLinearSolvers"
#include "Eigen/unsupported/Eigen/src/IterativeSolvers/GMRES.h"

// class HeatProblem;

class Inverser
{

public:
  Inverser(Parameters &p, HeatProblem &pb);
  Inverser(const Inverser &); // Declare copy constructor.
  Inverser &operator=(const Inverser &x);
  Inverser(){};
  ~Inverser();
  // Calculate the inverse algorithm
  double F(Eigen::VectorXd alpha,
           MatrixXd Immm,
           MatrixXd Ummm,
           double t,
           VectorXd sig,
           VectorXd sig0,
           VectorXd direction,
           VectorXd &angle_debut,
           VectorXd &angle_fin,
           VectorXd &Z);

  double F_TV(Eigen::VectorXd alpha,
              MatrixXd Immm,
              MatrixXd Ummm,
              double t,
              VectorXd sig,
              VectorXd sig0,
              VectorXd direction,
              VectorXd &angle_debut,
              VectorXd &angle_fin,
              VectorXd &Z);

  double cost_function(VectorXd &alpha,
                       VectorXd &angle_debut,
                       VectorXd &angle_fin,
                       VectorXd &sigma,
                       VectorXd &contact,
                       VectorXd &contact_terre,
                       MatrixXd &current_patern,
                       MatrixXd &measured_potentials,
                       double t,
                       VectorXd &dir_contact);

  double cost_function_both(VectorXd &alpha,
                            VectorXd &angle_debut,
                            VectorXd &angle_fin,
                            VectorXd &sigma,
                            VectorXd &sigma0,
                            VectorXd &contact,
                            VectorXd &contact_terre,
                            MatrixXd &current_patern,
                            MatrixXd &measured_potentials,
                            double t,
                            VectorXd &dir_contact,
                            VectorXd &dir_conductivity);

  void runinvcontact(VectorXd &alpha,
                     VectorXd &angle_debut,
                     VectorXd &angle_fin,
                     VectorXd &sigma);

  void runinv(VectorXd &alpha,
              VectorXd &angle_debut,
              VectorXd &angle_fin,
              VectorXd &Z);

  void runcombinedinv(VectorXd &alpha,
                      VectorXd &angle_debut,
                      VectorXd &angle_fin);

  double generator(float vect_meas);

private:
  // Initialization method
  // void startUp();

  Parameters *_params;
  Grid _grid;
  EITAssembler _assbl;
  Writer _writer;
  HeatProblem *_pb;
  Eigen::VectorXd _x, _y;
  double epsilon;

  // Pour trouver la direction.
  Eigen::VectorXd direction, sources, _sigma, _sigma0, _dir;
  double squaresize;

  SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> _solver;
  Eigen::GMRES<SparseMatrix<double>, Eigen::IncompleteLUT<double, int>> gmres;
};

#endif //_INVERSE_HPP_