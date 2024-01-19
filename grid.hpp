#ifndef _GRID_HPP_
#define _GRID_HPP_

#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"
#include "parameters.hpp"
#include <array>

// A cartesian grid. The size of the domain is fixed: it is [0,1]x[0,1]
class Grid
{

public:
  Grid(){};
  // Constructor with configuration (reads number of points in each direction)
  Grid(Parameters &p);

  // Returns the spacing discretization step in direction dim (0:x, 1:y)
  long double getSpacing(uint dim);

  // Returns the number of points in given dimension
  uint getNbPoints(uint dim);

  // return the radius of the geometry at hand
  double rayon(double theta, Eigen::VectorXd alpha);

  // the dirivative of the radius with respect to the angle
  long double rprime(double theta, Eigen::VectorXd alpha);
  // the dirivative of the radius with respect to the angle
  long double rdoubleprime(double theta, Eigen::VectorXd alpha);
  //
  double angle(double x, double y);

  //
  double rho(double theta, Eigen::VectorXd alpha);

  double f(double xx, double yy, double xc, double yc, Eigen::VectorXd alpha)
  {
    double r = std::sqrt(std::pow(xx - xc, 2) + std::pow(yy - yc, 2));
    double theta = std::atan2(yy - yc, xx - xc);
    return r - rayon(theta, alpha);
  }

  //
  double funcI(double theta, Eigen::VectorXd alpha);

  //
  double fintegral(double L, double theta_moin, double theta, Eigen::VectorXd alpha);

  //
  double angleplus(double L, double theta_moin, Eigen::VectorXd alpha);
  double angleplus_m(int i);

  // Gives the position in grid from matrix row index
  void index2D(uint row, uint &i, uint &j);

  // Get coordinates
  void coords(uint i, uint j, double &x, double &y);

  // Get coordinates
  long double Levelset(uint i, uint j, Eigen::VectorXd alpha);

  // Construit le vecteur de numerotation des point irregulier
  void BuildT(Eigen::VectorXd alpha);

  // Gives the matrix row ID given position in grid  in irregular counting
  uint ind(uint i, uint j);

  // Dertrmnination du numéro de l'electrode auquel appartient le point (x_i,y_j)
  double Funczone(double thetha, Eigen::VectorXd alpha, Eigen::VectorXd angle_debut, Eigen::VectorXd angle_fin);

  // Construit la derivé de la level set
  void DLevelset(Eigen::VectorXd alpha);

  // Construit la fonction Delta pour discretisé l'integrale
  double Delta(int i, int j, Eigen::VectorXd alpha);

  // retrouve les point d'interface
  void inter(Eigen::VectorXd alpha);

  // retrouve les point d'interface
  void Zonelectrode(Eigen::VectorXd alpha, Eigen::VectorXd angle_debut, Eigen::VectorXd angle_fin);

  // conductivity function
  double coef(double x, double y);

  //
  double sign_fonc(double x);

  //
  void discret_scheme(uint i, uint j, double fix, double fiy, double x, double y, Eigen::MatrixXd &M, Eigen::MatrixXd &M_ind);
  //
  long double source(double x, double y);

  //
  long double val_du(double x, double y, double fix, double fiy);
  //
  long double val_dt(double x, double y, Eigen::VectorXd alpha);
  //
  long double val_dt_norm(double x, double y, Eigen::VectorXd alpha);
  //
  long double cond_inter(double x, double y);
  //
  long double normalvec_x(double thetha, Eigen::VectorXd alpha);
  //
  long double normalvec_y(double thetha, Eigen::VectorXd alpha);

  long double tangentvec_x(double thetha, Eigen::VectorXd alpha);
  //
  long double tangentvec_y(double thetha, Eigen::VectorXd alpha);
  //
  double angle_moin(int i);
  //
  double courbure(double thetha, Eigen::VectorXd alpha);
  //
  void ordre();
  //
  void coord_int();
  //
  double longeur_arc(double theta1, double theta2, Eigen::VectorXd alpha);
  //
  long double tangent_vec(double x, double y, Eigen::VectorXd alpha);

  // Pour récuperer la matrice du laplacien
  int get_tx2(int i, int j) { return _tx2(i, j); };
  int get_ty2(int i, int j) { return _ty2(i, j); };
  int get_xix2_2(int i, int j) { return _xix2_2(i, j); };
  int get_xiy2_2(int i, int j) { return _xiy2_2(i, j); };
  double get_xinter(int i, int j) { return xinter(i, j); };
  double get_yinter(int i, int j) { return yinter(i, j); };
  double get_u(int i, int j) { return _u(i, j); };
  double get_v(int i, int j) { return _v(i, j); };
  double get_xi(int i, int j) { return x_i; };
  double get_yi(int i, int j) { return y_i; };
  double get_xk(int i, int j) { return x_k; };
  double get_yk(int i, int j) { return y_k; };
  int get_indk(int i, int j) { return indice_k; };
  int get_indi(int i, int j) { return indice_i; };
  double get_x0(int k) { return x_0(k); };
  double get_y0(int k) { return y_0(k); };
  double get_x2(int k) { return x_2(k); };
  double get_y2(int k) { return y_2(k); };
  int get_ind0(int k) { return ind_0(k); };
  int get_ind2(int k) { return ind_2(k); };
  double get_angle_ind(int k) { return angles_interface(0, k); };
  double get_angle_ind_x(int k) { return angles_interface(1, k); };
  double get_angle_ind_y(int k) { return angles_interface(2, k); };
  int get_inter_ind(int k) { return indice_interface(k); };
  int get_inter_fin(int k) { return indice_final(k); };
  double get_coord_int(int i, int k) { return coord_inter(i, k); };

  Eigen::VectorXd Get_zone() { return zone; };
  int get_zone(int i) { return zone(i); };
  double Get_nbsigne() { return _nbsigne; };

private:
  std::array<uint, 3> _N;
  std::array<double, 2> _h;

  int _r, _nbsigne;
  double _dxphiplus, _dxphimoins, _dxphizero;
  double _dyphiplus, _dyphimoins, _dyphizero;
  double _absgrad;
  double _delta;
  double _alpha;

  // variable pour determiner le shema de discretisa

  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> _ty2, _tx2;
  Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> _xiy2_2, _xix2_2;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _u, _v;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> xinter, yinter;
  Eigen::VectorXd _x, _y;
  Eigen::VectorXd zone;
  double x_k, y_k, x_i, y_i;
  int indice_k, indice_i;
  Eigen::VectorXd x_0, y_0, x_2, y_2;
  Eigen::VectorXi ind_0, ind_2;
  Eigen::MatrixXd coord_inter;
  Eigen::MatrixXd angles_interface;
  Eigen::VectorXi indice_interface, indice_final;
  double squaresize;
};

#endif //_GRID_HPP_
