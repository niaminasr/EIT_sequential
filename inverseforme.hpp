#ifndef _INVERSEFORME_HPP_
#define _INVERSEFORME_HPP_


#include "parameters.hpp"
#include "grid.hpp"
#include "assembler.hpp"
#include "assembler_forme.hpp"
#include "writer.hpp"
#include "heatProblem.hpp"
#include "inverse.hpp"
#include "Eigen/Eigen/IterativeLinearSolvers"
#include "Eigen/unsupported/Eigen/src/IterativeSolvers/GMRES.h"


class Inverserr {

  public:

    // Initializes the problem from a text input file
    Inverserr(Parameters& p, HeatProblem& pb);
    //runing the inverse problem
    void run_inv_cond_in_any_geom();

    void run_inv_form();

    void run_inv_electrode_position();
  
  private:

    // Initialization method
    Parameters*   _params;
    Grid          _grid;
    EITAssembler _assbl;
    HeatAssemblerforme _assblf;
    Writer        _writer;
    Inverser* _inv;
    HeatProblem*  _pb;
 
    Eigen::VectorXd _x, _y;
    Eigen::VectorXd _solf;
    Eigen::MatrixXd _R,_b,_RR,_Rm;
    double dx,dy;
    double epsilon;

    double squaresize;
  
  

    //Pour trouver la direction.
    double direction, F;
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > _solver;
    Eigen::MatrixXd _dxsoll,_dysoll;
    double exactsol,linf,errmax,resmax,xp,xm,uup,uum,errdxmax,dxsolexact;
};

#endif//_INVERSEFORM_HPP_
