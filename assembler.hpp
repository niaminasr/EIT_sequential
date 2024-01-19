#ifndef ASSEMBLER
#define ASSEMBLER

#include "grid.hpp"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Eigen"
#include "Eigen/Eigen/Sparse"
#include <iostream>

using namespace Eigen;

class HeatProblem;

class EITAssembler
{

public:
  EITAssembler(){};
  EITAssembler(Grid *g, HeatProblem *);

  /*Creates the lines related to the jump conditions of the
  interface points on the interface: either we are on an electrode or not*/
  void fluxx(int i, int j);
  void fluxy(int i, int j);
  void electrodefluxx(int i, int j, VectorXd Z);
  void electrodefluxy(int i, int j, VectorXd Z);

  // Build the Matrix of the EIT problem :: we assemble
  void BuildLaplacianMatrix(VectorXd alpha,
                            VectorXd angle_debut,
                            VectorXd angle_fin,
                            VectorXd sigma,
                            VectorXd Z);

  // Builds the EIT source term
  void BuildSourceTerm(VectorXd &alpha,
                       VectorXd &Imm,
                       VectorXd &Z);
  // Builds the manufactured source term
  void BuildSourceTermMAN(VectorXd &alpha, VectorXd &Z);
  // builds the exact solution
  void assemblesolexact();
  // Clearing function
  void CLEAR();

  void BuildDirMatrix(VectorXd alpha,
                      VectorXd angle_debut,
                      VectorXd angle_fin,
                      VectorXd sigma);
  // buils the first part of sourse term
  void Buildsourcedir1(uint l, MatrixXd u, MatrixXd w, VectorXd& alpha);
  // builds the second part of the source term
  void Buildsourcedir2(Eigen::VectorXd u, Eigen::VectorXd v, VectorXd alpha);

  VectorXd grad_of_vec(VectorXd alpha, VectorXd vec);

  VectorXd Build_dirVT_G(VectorXd alpha, VectorXd vec);
  // normal derivative
  double grad_normal_x(Eigen::MatrixXd &u, int i, int j, int l);
  // normal derivative
  double grad_normal_y(Eigen::MatrixXd &u, int i, int j, int l);

  // Get functions ::
  Eigen::SparseMatrix<double> GetMatrix() { return _lap_mat; };
  Eigen::SparseMatrix<double> GetMatrixdir() { return _lap_dir; };
  Eigen::VectorXd GetSourceTerm() { return _source_term; };
  Eigen::VectorXd GetSourceTermMan() { return _source_term_man; };
  double Getdsigma(uint i, uint j) { return dsigma(i, j); };
  double Getlapdir(uint i, uint j) { return lap_dir(i, j); };


  double Compute_integral(MatrixXd &u,
                          MatrixXd &w,
                          VectorXd &alpha,
                          VectorXd &angle_debut,
                          VectorXd &angle_fin,
                          int i,
                          int l);

private:
  HeatProblem *_pb; // Source term function
  Grid *_grid;

  Eigen::SparseMatrix<double> _lap_mat, _lap_dir;
  Eigen::VectorXd _source_term, _source_term_man;

  typedef Triplet<double> Trip;
  std::vector<Trip> trp;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dsigma, lap_dir;

  VectorXd _x, _y;
  VectorXd mat3, matt3;
  VectorXd mat2, matt2;
  VectorXd mat1, matt1;

  VectorXd Im;
  VectorXd Um;
  VectorXd Ummm;

  Eigen::VectorXd sigma;
  int m_size;
  int Nx, Ny, Ne, Ni, n;
  double dx, dy;
  int marqueur, indk, indi;
  double rr, fix, fiy, theta, fiix, fiiy, d, dd;
  double alphaj, betaj;
  double alphak, betak;
  double alphai, betai, denom;
  double xi, yi, xk, yk, xj, yj;
  double pi;
  int m;

  double squaresize;
};

#endif /* ASSEMBLER */
