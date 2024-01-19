#ifndef _HEAT_PROBLEM_HPP_
#define _HEAT_PROBLEM_HPP_

#include "parameters.hpp"
#include "grid.hpp"
#include "assembler.hpp"
#include "assembler_forme.hpp"
#include "writer.hpp"
#include "Eigen/Eigen/IterativeLinearSolvers"
#include "Eigen/unsupported/Eigen/src/IterativeSolvers/GMRES.h"

// A class that solve the linear system on a cartesian grid
class HeatProblem
{

public:
  // Initializes the problem from a text input file
  HeatProblem(Parameters &p);
  // Does the whole computation of the solution
  void run(VectorXd _sigma,
           VectorXd &alpha,
           VectorXd &Imm,
           VectorXd &angle_debut,
           VectorXd &angle_fin,
           VectorXd &Z);
  //
  void runwithsource(VectorXd &_sigma,
                     VectorXd &alpha,
                     VectorXd &angle_debut,
                     VectorXd &angle_fin,
                     VectorXd &Z);
  // source
  double sourceTerm();
  // Diffusion coefficient is scalar in this case
  double getDiffusionCoeff() { return _k; };
  // Récupère la Solution
  Eigen::VectorXd GetSolution() { return _sol; };

private:
  // Initialization method
  // void startUp();

  Parameters *_params;
  Grid _grid;
  EITAssembler _assbl;
  Writer _writer;
  Eigen::VectorXd _x, _y, _alpha;

  double _k;
  // Solution
  Eigen::VectorXd _sol, residu, yy, xx, xxx;
  Eigen::VectorXd sigma;
  // Solver
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> _solver;
  Eigen::GMRES<SparseMatrix<double>, Eigen::IncompleteLUT<double, int>> gmres;
  Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double, int>> bicg;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _dxsol, _dysol;

  double exactsol, linf, errmax, resmax, xp, xm, uup, uum, errdxmax, dxsolexact;
  double squaresize;
};

#endif //_HEAT_PROBLEM_HPP_
