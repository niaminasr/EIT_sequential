#ifndef _ASSEMBLER_FORME_HPP_
#define _ASSEMBLER_FORME_HPP_

#include "grid.hpp"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Eigen"
#include "Eigen/Eigen/Sparse"
#include <iostream>

using namespace Eigen;

class HeatProblem;

class HeatAssemblerforme{

  public :

    HeatAssemblerforme(){};
    HeatAssemblerforme(Grid* g, HeatProblem*);

    //Fills the matrix if the interface point is  not on an electode
    void fluxx_forme(int i, int j);
    //Fills the matrix if the interface point is  not on an electode
    void fluxy_forme(int i, int j);
    //Fills the matrix if the interface point is on an electrode
    void electrodefluxx_forme(int i, int j);
    //Fills the matrix if the interface point is on an electrode
    void electrodefluxy_forme(int i, int j);
    //assembles laplacien Matrix
    void BuildLaplacianMatrix_forme(Eigen::VectorXd alpha);
    // Construit le terme source EIT
    void BuildSourceTerm_forme(Eigen::VectorXd& alpha,Eigen::MatrixXd& u, Eigen::VectorXd& Imm);
	  void BuildSourceTerm_forme_source(Eigen::VectorXd& alpha, Eigen::MatrixXd& u);
    //builds the exact solution
    void assemblesolexact();
    //Clearing function
    void CLEAR();
    
    // Pour récuperer la matrice du laplacien
	Eigen::SparseMatrix<double> GetMatrixf(){return _lap_mat;};
	// Pour récupérer le terme source
	Eigen::VectorXd  GetSourceTermf(){return _source_term;};
    // Pour récupérer le terme source
	Eigen::VectorXd  GetSourceTermfn(){return _source_term2;};



  private:

    HeatProblem* _pb; // Source term function
    Grid*        _grid;


    Eigen::SparseMatrix<double> _lap_mat;
  	Eigen::VectorXd _source_term,_source_term2;

    typedef Triplet<double> Trip;
    std::vector<Trip> trp;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> lap_dir;
    Eigen::MatrixXf a,b;
    //vector<Triplet<double>> triplets;
  	// vecteur x = {xmin + hx, ... , xmax - hx}
  	// vecteur y = {ymin + hy, ... , ymax - hy}
  	Eigen::VectorXd _x, _y;
    //const StorageIndex& Eigen::Triplet<Scalar, StorageIndex>::col()




    Eigen::VectorXd mat3, matt3;
    Eigen::VectorXd mat2, matt2;
    Eigen::VectorXd mat1, matt1,qqq;
    Eigen::VectorXd Im;
    Eigen::VectorXd Um;
    Eigen::VectorXd Ummm;

    Eigen::VectorXd sigma;
    int Nx,Ny,Ne,n;
    double dx,dy;
    int marqueur,indk,indi;
    double rr, fix, fiy, theta,fiix,fiiy,d,dd;
    double alphak, betak, xi,yi,xk,yk,xj,yj;
    double alphai, betai ,denom ,pi;
    double alphaj, betaj ;
    double r2,the,_eta;
    int m;
    int r_0,somme,dthetaa,dth,Lon;
    double coeff;

    double squaresize;

};

#endif // _ASSEMBLERFORME_HPP_
