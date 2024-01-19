#include "assembler.hpp"
#include "heatProblem.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "Eigen/Eigen/Eigen"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

using namespace std;
using namespace Eigen;

EITAssembler::EITAssembler(Grid *g, HeatProblem *pb)
{

  _grid = g;
  _pb = pb;

  Nx = _grid->getNbPoints(0);
  Ny = _grid->getNbPoints(1);
  Ne = _grid->getNbPoints(2);
  Ni = _grid->Get_nbsigne();

  dx = _grid->getSpacing(0);
  dy = _grid->getSpacing(1);

  mat1.resize(8 * Nx * Ny);
  mat2.resize(8 * Nx * Ny);
  mat3.resize(8 * Nx * Ny);
  Im.resize(Ne);
  Um.resize(Ne);
  pi = acos(-1.);

  squaresize = 2.;

  _x.resize(Nx); // (xmin+h, ..., xmax-h)
  for (int i = 0; i < Nx; ++i)
    _x[i] = ((i + 1) - 0.5) * dx - squaresize;

  _y.resize(Ny);
  for (int i = 0; i < Ny; ++i)
    _y[i] = ((i + 1) - 0.5) * dy - squaresize;
}

void EITAssembler::fluxx(int i, int j)
{

  fix = 0;
  fiy = 0;

  // calcul normale sur le point d'interface
  fix = (_grid->get_u(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_u(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
  fiy = (_grid->get_v(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_v(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
  fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
  fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
  fix = fiix;
  fiy = fiiy;

  Eigen::MatrixXd M, M_ind;
  M.resize(6, 2);
  M_ind.resize(6, 2);

  M(0, 0) = _x[i];
  M(0, 1) = _y[j + 1];
  M(1, 0) = _x[i + 1];
  M(1, 1) = _y[j + 1];
  M(2, 0) = _x[i];
  M(2, 1) = _y[j];
  M(3, 0) = _x[i + 1];
  M(3, 1) = _y[j];
  M(4, 0) = _x[i];
  M(4, 1) = _y[j - 1];
  M(5, 0) = _x[i + 1];
  M(5, 1) = _y[j - 1];

  M_ind(0, 0) = i;
  M_ind(0, 1) = j + 1;
  M_ind(1, 0) = i + 1;
  M_ind(1, 1) = j + 1;
  M_ind(2, 0) = i;
  M_ind(2, 1) = j;
  M_ind(3, 0) = i + 1;
  M_ind(3, 1) = j;
  M_ind(4, 0) = i;
  M_ind(4, 1) = j - 1;
  M_ind(5, 0) = i + 1;
  M_ind(5, 1) = j - 1;

  xj = _grid->get_xinter(i, j);
  yj = _y[j];
  _grid->discret_scheme(i, j, fix, fiy, xj, yj, M, M_ind);

  xk = _grid->get_xk(i, j);
  yk = _grid->get_yk(i, j);
  indk = _grid->get_indk(i, j);

  xi = _grid->get_xi(i, j);
  yi = _grid->get_yi(i, j);
  indi = _grid->get_indi(i, j);

  // calcul des coef du gradient (interpolation lineaire sur les 3 points)
  denom = (xj - xk) * (yj - yi) - (xj - xi) * (yj - yk);
  alphaj = (yk - yi) / denom;
  betaj = (xi - xk) / denom;

  denom = (xi - xk) * (yi - yj) - (xi - xj) * (yi - yk);
  alphai = (yk - yj) / denom;
  betai = (xj - xk) / denom;

  denom = (xk - xj) * (yk - yi) - (xk - xi) * (yk - yj);
  alphak = (yj - yi) / denom;
  betak = (xi - xj) / denom;

  // indk
  mat1[m] = _grid->get_tx2(i, j) - 1;
  mat2[m] = indk - 1;
  mat3[m] = _grid->coef(_grid->get_xinter(i, j), _y[j]) * (alphak * fix + betak * fiy);

  m = m + 1;

  // indi
  mat1[m] = _grid->get_tx2(i, j) - 1;
  mat2[m] = indi - 1;
  mat3[m] = _grid->coef(_grid->get_xinter(i, j), _y[j]) * (alphai * fix + betai * fiy);
  m = m + 1;

  // int(i,i+1),j
  mat1[m] = _grid->get_tx2(i, j) - 1;
  mat2[m] = _grid->get_tx2(i, j) - 1;
  mat3[m] = _grid->coef(_grid->get_xinter(i, j), _y[j]) * (alphaj * fix + betaj * fiy);
  m = m + 1;
}

void EITAssembler::fluxy(int i, int j)
{

  fix = 0;
  fiy = 0;

  // calcul normale sur le point d'interface
  fix = (_grid->get_u(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_u(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
  fiy = (_grid->get_v(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_v(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
  fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
  fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
  fix = fiix;
  fiy = fiiy;

  Eigen::MatrixXd M, M_ind;
  M.resize(6, 2);
  M_ind.resize(6, 2);

  M(0, 0) = _x[i - 1];
  M(0, 1) = _y[j + 1];
  M(1, 0) = _x[i - 1];
  M(1, 1) = _y[j];
  M(2, 0) = _x[i];
  M(2, 1) = _y[j + 1];
  M(3, 0) = _x[i];
  M(3, 1) = _y[j];
  M(4, 0) = _x[i + 1];
  M(4, 1) = _y[j + 1];
  M(5, 0) = _x[i + 1];
  M(5, 1) = _y[j];

  M_ind(0, 0) = i - 1;
  M_ind(0, 1) = j + 1;
  M_ind(1, 0) = i - 1;
  M_ind(1, 1) = j;
  M_ind(2, 0) = i;
  M_ind(2, 1) = j + 1;
  M_ind(3, 0) = i;
  M_ind(3, 1) = j;
  M_ind(4, 0) = i + 1;
  M_ind(4, 1) = j + 1;
  M_ind(5, 0) = i + 1;
  M_ind(5, 1) = j;

  xj = _x[i];
  yj = _grid->get_yinter(i, j);
  _grid->discret_scheme(i, j, fix, fiy, xj, yj, M, M_ind);

  xk = _grid->get_xk(i, j);
  yk = _grid->get_yk(i, j);
  indk = _grid->get_indk(i, j);

  xi = _grid->get_xi(i, j);
  yi = _grid->get_yi(i, j);
  indi = _grid->get_indi(i, j);

  // calcul des coef du gradient (interpolation lineaire sur les 3 points)
  denom = (xj - xk) * (yj - yi) - (xj - xi) * (yj - yk);
  alphaj = (yk - yi) / denom;
  betaj = (xi - xk) / denom;

  denom = (xi - xk) * (yi - yj) - (xi - xj) * (yi - yk);
  alphai = (yk - yj) / denom;
  betai = (xj - xk) / denom;

  denom = (xk - xj) * (yk - yi) - (xk - xi) * (yk - yj);
  alphak = (yj - yi) / denom;
  betak = (xi - xj) / denom;

  // indk
  mat1[m] = _grid->get_ty2(i, j) - 1;
  mat2[m] = indk - 1;
  mat3[m] = _grid->coef(_x[i], _grid->get_yinter(i, j)) * (alphak * fix + betak * fiy);
  m = m + 1;

  // indi
  mat1[m] = _grid->get_ty2(i, j) - 1;
  mat2[m] = indi - 1;
  mat3[m] = _grid->coef(_x[i], _grid->get_yinter(i, j)) * (alphai * fix + betai * fiy);
  m = m + 1;

  // int(i,i+1),j
  mat1[m] = _grid->get_ty2(i, j) - 1;
  mat2[m] = _grid->get_ty2(i, j) - 1;
  mat3[m] = _grid->coef(_x[i], _grid->get_yinter(i, j)) * (alphaj * fix + betaj * fiy);
  m = m + 1;
}

void EITAssembler::electrodefluxx(int i, int j, VectorXd Z)
{

  fix = 0;
  fiy = 0;

  // std::cout << Z.size() << std::endl;

  // calcul normale sur le point d'interface
  fix = (_grid->get_u(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_u(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
  fiy = (_grid->get_v(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_v(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
  fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
  fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
  fix = fiix;
  fiy = fiiy;

  Eigen::MatrixXd M, M_ind;
  M.resize(6, 2);
  M_ind.resize(6, 2);

  M(0, 0) = _x[i];
  M(0, 1) = _y[j + 1];
  M(1, 0) = _x[i + 1];
  M(1, 1) = _y[j + 1];
  M(2, 0) = _x[i];
  M(2, 1) = _y[j];
  M(3, 0) = _x[i + 1];
  M(3, 1) = _y[j];
  M(4, 0) = _x[i];
  M(4, 1) = _y[j - 1];
  M(5, 0) = _x[i + 1];
  M(5, 1) = _y[j - 1];

  M_ind(0, 0) = i;
  M_ind(0, 1) = j + 1;
  M_ind(1, 0) = i + 1;
  M_ind(1, 1) = j + 1;
  M_ind(2, 0) = i;
  M_ind(2, 1) = j;
  M_ind(3, 0) = i + 1;
  M_ind(3, 1) = j;
  M_ind(4, 0) = i;
  M_ind(4, 1) = j - 1;
  M_ind(5, 0) = i + 1;
  M_ind(5, 1) = j - 1;

  xj = _grid->get_xinter(i, j);
  yj = _y[j];
  _grid->discret_scheme(i, j, fix, fiy, xj, yj, M, M_ind);

  xk = _grid->get_xk(i, j);
  yk = _grid->get_yk(i, j);
  indk = _grid->get_indk(i, j);

  xi = _grid->get_xi(i, j);
  yi = _grid->get_yi(i, j);
  indi = _grid->get_indi(i, j);

  // calcul des coef du gradient (interpolation lineaire sur les 3 points)
  denom = (xj - xk) * (yj - yi) - (xj - xi) * (yj - yk);
  alphaj = (yk - yi) / denom;
  betaj = (xi - xk) / denom;

  denom = (xi - xk) * (yi - yj) - (xi - xj) * (yi - yk);
  alphai = (yk - yj) / denom;
  betai = (xj - xk) / denom;

  denom = (xk - xj) * (yk - yi) - (xk - xi) * (yk - yj);
  alphak = (yj - yi) / denom;
  betak = (xi - xj) / denom;

  // indk
  mat1[m] = _grid->get_tx2(i, j) - 1;
  mat2[m] = indk - 1;
  mat3[m] = _grid->coef(_grid->get_xinter(i, j), _y[j]) * (alphak * fix + betak * fiy);
  m = m + 1;

  // indi
  mat1[m] = _grid->get_tx2(i, j) - 1;
  mat2[m] = indi - 1;
  mat3[m] = _grid->coef(_grid->get_xinter(i, j), _y[j]) * (alphai * fix + betai * fiy);
  m = m + 1;

  // indj
  mat1[m] = _grid->get_tx2(i, j) - 1;
  mat2[m] = _grid->get_tx2(i, j) - 1;
  mat3[m] = Z(_grid->get_zone(_grid->get_tx2(i, j)) - 1) * 1 + _grid->coef(_grid->get_xinter(i, j), _y[j]) * (alphaj * fix + betaj * fiy);
  m = m + 1;

  //
  mat1[m] = _grid->get_tx2(i, j) - 1;
  mat2[m] = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i, j)) - 1;
  mat3[m] = -1. * Z(_grid->get_zone(_grid->get_tx2(i, j)) - 1);
  m = m + 1;
}

void EITAssembler::electrodefluxy(int i, int j, VectorXd Z)
{

  int l, k;
  fix = 0;
  fiy = 0;

  // calcul normale sur le point d'interface
  fix = (_grid->get_u(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_u(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
  fiy = (_grid->get_v(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_v(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
  fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
  fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
  fix = fiix;
  fiy = fiiy;

  Eigen::MatrixXd M, M_ind;
  M.resize(6, 2);
  M_ind.resize(6, 2);

  M(0, 0) = _x[i - 1];
  M(0, 1) = _y[j + 1];
  M(1, 0) = _x[i - 1];
  M(1, 1) = _y[j];
  M(2, 0) = _x[i];
  M(2, 1) = _y[j + 1];
  M(3, 0) = _x[i];
  M(3, 1) = _y[j];
  M(4, 0) = _x[i + 1];
  M(4, 1) = _y[j + 1];
  M(5, 0) = _x[i + 1];
  M(5, 1) = _y[j];

  M_ind(0, 0) = i - 1;
  M_ind(0, 1) = j + 1;
  M_ind(1, 0) = i - 1;
  M_ind(1, 1) = j;
  M_ind(2, 0) = i;
  M_ind(2, 1) = j + 1;
  M_ind(3, 0) = i;
  M_ind(3, 1) = j;
  M_ind(4, 0) = i + 1;
  M_ind(4, 1) = j + 1;
  M_ind(5, 0) = i + 1;
  M_ind(5, 1) = j;

  xj = _x[i];
  yj = _grid->get_yinter(i, j);
  _grid->discret_scheme(i, j, fix, fiy, xj, yj, M, M_ind);

  xk = _grid->get_xk(i, j);
  yk = _grid->get_yk(i, j);
  indk = _grid->get_indk(i, j);

  xi = _grid->get_xi(i, j);
  yi = _grid->get_yi(i, j);
  indi = _grid->get_indi(i, j);

  // calcul des coef du gradient (interpolation lineaire sur les 3 points)
  denom = (xj - xk) * (yj - yi) - (xj - xi) * (yj - yk);
  alphaj = (yk - yi) / denom;
  betaj = (xi - xk) / denom;

  denom = (xi - xk) * (yi - yj) - (xi - xj) * (yi - yk);
  alphai = (yk - yj) / denom;
  betai = (xj - xk) / denom;

  denom = (xk - xj) * (yk - yi) - (xk - xi) * (yk - yj);
  alphak = (yj - yi) / denom;
  betak = (xi - xj) / denom;

  // indk
  mat1[m] = _grid->get_ty2(i, j) - 1;
  mat2[m] = indk - 1;
  mat3[m] = _grid->coef(_x[i], _grid->get_yinter(i, j)) * (alphak * fix + betak * fiy);
  m = m + 1;

  // indi
  mat1[m] = _grid->get_ty2(i, j) - 1;
  mat2[m] = indi - 1;
  mat3[m] = _grid->coef(_x[i], _grid->get_yinter(i, j)) * (alphai * fix + betai * fiy);
  m = m + 1;

  // int(i,i+1),j
  mat1[m] = _grid->get_ty2(i, j) - 1;
  mat2[m] = _grid->get_ty2(i, j) - 1;
  mat3[m] = Z(_grid->get_zone(_grid->get_ty2(i, j)) - 1) * 1 + _grid->coef(_x[i], _grid->get_yinter(i, j)) * (alphaj * fix + betaj * fiy);
  m = m + 1;

  // Um
  mat1[m] = _grid->get_ty2(i, j) - 1;
  mat2[m] = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i, j)) - 1;
  mat3[m] = -1. * Z(_grid->get_zone(_grid->get_ty2(i, j)) - 1);
  m = m + 1;
}

void EITAssembler::BuildLaplacianMatrix(VectorXd alpha,
                                        VectorXd angle_debut,
                                        VectorXd angle_fin,
                                        VectorXd sigma,
                                        VectorXd Z)
{

  double b;

  _grid->BuildT(alpha);
  _grid->inter(alpha);
  _grid->DLevelset(alpha);
  _grid->Zonelectrode(alpha, angle_debut, angle_fin);
  int size = Nx * Ny + _grid->Get_nbsigne() + Ne;

  _lap_mat.resize(size, size);
  sigma.resize(size);

  double coeff = 0;
  m = 0;

  // Premier cas:: i=0; j=0.
  mat1[m] = _grid->ind(0, 0) - 1;
  mat2[m] = _grid->ind(0, 0) - 1;
  mat3[m] = 0;

  coeff = _grid->coef(_x[0] + dx / 2, _y[0]);
  if (_grid->get_xix2_2(0, 0) != 0)
    coeff = _grid->coef(_grid->get_xinter(0, 0) / 2 + _x[0] / 2, _y[0]);
  mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(0, 0)) + _grid->get_xix2_2(0, 0) * abs(_grid->get_xinter(0, 0) - _x[0])) / dx);

  coeff = _grid->coef(_x[0], _y[0]) / 2;
  mat3[m] = mat3[m] + coeff * (-1. / (dx * (1.)) / dx);

  coeff = _grid->coef(_x[0], _y[0] + dy / 2.);
  if (_grid->get_xiy2_2(0, 0) != 0)
    coeff = _grid->coef(_x[0], _grid->get_yinter(0, 0) / 2 + _y[0] / 2);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, 0)) + _grid->get_xiy2_2(0, 0) * abs(_grid->get_yinter(0, 0) - _y[0])) / dx);

  coeff = _grid->coef(_x[0], _y[0]) / 2;
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1.)) / dx);
  m = m + 1;

  // i=0, j=Ny-1
  mat1[m] = _grid->ind(0, Ny - 1) - 1;
  mat2[m] = _grid->ind(0, Ny - 1) - 1;
  mat3[m] = 0;

  coeff = _grid->coef(_x[0] + dx / 2, _y[Ny - 1]);
  if (_grid->get_xix2_2(0, Ny - 1) != 0)
    coeff = _grid->coef(_grid->get_xinter(0, Ny - 1) / 2 + _x[0] / 2, _y[Ny - 1]);
  mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(0, Ny - 1)) + _grid->get_xix2_2(0, Ny - 1) * abs(_grid->get_xinter(0, Ny - 1) - _x[0])) / dx);

  coeff = _grid->coef(_x[0] - dx / 2, _y[Ny - 1]);
  mat3[m] = mat3[m] + coeff * (-1. / (dx * (1.)) / dx);

  coeff = _grid->coef(_x[0], _y[Ny - 1]) / 2;
  if (_grid->get_xiy2_2(0, Ny - 1) != 0)
    coeff = _grid->coef(_x[0], _grid->get_yinter(0, Ny - 1) / 2 + _y[Ny - 1] / 2);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, Ny - 1)) + _grid->get_xiy2_2(0, Ny - 1) * abs(_grid->get_yinter(0, Ny - 1) - _y[Ny - 1])) / dx);

  coeff = _grid->coef(_x[0], _y[Ny - 1] - dy / 2);
  if (_grid->get_xiy2_2(0, Ny - 1 - 1) != 0)
    coeff = _grid->coef(_x[0], _grid->get_yinter(0, Ny - 1 - 1) / 2. + _y[Ny - 1] / 2.);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, Ny - 1 - 1)) + _grid->get_xiy2_2(0, Ny - 1 - 1) * abs(_y[Ny - 1] - _grid->get_yinter(0, Ny - 1 - 1))) / dx);
  m = m + 1;

  // i=Nx-1, j=0
  mat1[m] = _grid->ind(Nx - 1, 0) - 1;
  mat2[m] = _grid->ind(Nx - 1, 0) - 1;
  mat3[m] = 0;

  coeff = _grid->coef(_x[Nx - 1] + dx / 2, _y[0]);
  if (_grid->get_xix2_2(Nx - 1, 0) != 0)
    coeff = _grid->coef(_grid->get_xinter(Nx - 1, 0) / 2 + _x[Nx - 1] / 2, _y[0]);
  mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(Nx - 1, 0)) + _grid->get_xix2_2(Nx - 1, 0) * abs(_grid->get_xinter(Nx - 1, 0) - _x[Nx - 1])) / dx);

  coeff = _grid->coef(_x[Nx - 1] - dx / 2, _y[0]);
  if (_grid->get_xix2_2(Nx - 1 - 1, 0) != 0)
    coeff = _grid->coef(_grid->get_xinter(Nx - 1 - 1, 0) / 2 + _x[Nx - 1] / 2, _y[0]);
  mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(Nx - 1 - 1, 0)) + _grid->get_xix2_2(Nx - 1 - 1, 0) * abs(_x[Nx - 1] - _grid->get_xinter(Nx - 1 - 1, 0))) / dx);

  coeff = _grid->coef(_x[Nx - 1], _y[0] + dy / 2.);
  if (_grid->get_xiy2_2(Nx - 1, 0) != 0)
    coeff = _grid->coef(_x[Nx - 1], _grid->get_yinter(Nx - 1, 0) / 2 + _y[0] / 2);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(Nx - 1, 0)) + _grid->get_xiy2_2(Nx - 1, 0) * abs(_grid->get_yinter(Nx - 1, 0) - _y[0])) / dx);

  coeff = _grid->coef(_x[Nx - 1], _y[0]) / 2;
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1.)) / dx);
  m = m + 1;

  // i=Nx-1, j=Ny-1
  mat1[m] = _grid->ind(Nx - 1, Ny - 1) - 1;
  mat2[m] = _grid->ind(Nx - 1, Ny - 1) - 1;
  mat3[m] = 0;

  int i = Nx - 1;
  int j = Ny - 1;

  coeff = _grid->coef(_x[i], _y[j]) / 2;
  if (_grid->get_xix2_2(i, j) != 0)
    coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
  mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx);

  coeff = _grid->coef(_x[i] - dx / 2, _y[j]);
  if (_grid->get_xix2_2(i - 1, j) != 0)
    coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2 + _x[i] / 2, _y[j]);
  mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i - 1, j)) + _grid->get_xix2_2(i - 1, j) * abs(_x[i] - _grid->get_xinter(i - 1, j))) / dx);

  coeff = _grid->coef(_x[i], _y[j]) / 2;
  if (_grid->get_xiy2_2(i, j) != 0)
    coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

  coeff = _grid->coef(_x[i], _y[j] - dy / 2);
  if (_grid->get_xiy2_2(i, j - 1) != 0)
    coeff = _grid->coef(_x[i], _grid->get_yinter(i, j - 1) / 2. + _y[j] / 2.);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j - 1)) + _grid->get_xiy2_2(i, j - 1) * abs(_y[j] - _grid->get_yinter(i, j - 1))) / dx);
  m = m + 1;

  for (int i = 1; i < Nx - 1; ++i)
  {
    int j = 0;
    mat1[m] = _grid->ind(i, j) - 1;
    mat2[m] = _grid->ind(i, j) - 1;
    mat3[m] = 0;

    coeff = _grid->coef(_x[i] + dx / 2, _y[j]);
    if (_grid->get_xix2_2(i, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx);

    coeff = _grid->coef(_x[i] - dx / 2, _y[j]);
    if (_grid->get_xix2_2(i - 1, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i - 1, j)) + _grid->get_xix2_2(i - 1, j) * abs(_x[i] - _grid->get_xinter(i - 1, j))) / dx);

    coeff = _grid->coef(_x[i], _y[j] + dy / 2.);
    if (_grid->get_xiy2_2(i, j) != 0)
      coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

    coeff = _grid->coef(_x[i], _y[j]) / 2;
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1.)) / dx);
    m = m + 1;
  }

  for (int i = 1; i < Nx - 1; ++i)
  {
    int j = Ny - 1;
    mat1[m] = _grid->ind(i, j) - 1;
    mat2[m] = _grid->ind(i, j) - 1;
    mat3[m] = 0;

    coeff = _grid->coef(_x[i] + dx / 2, _y[j]);
    if (_grid->get_xix2_2(i, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx);

    coeff = _grid->coef(_x[i] - dx / 2, _y[j]);
    if (_grid->get_xix2_2(i - 1, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i - 1, j)) + _grid->get_xix2_2(i - 1, j) * abs(_x[i] - _grid->get_xinter(i - 1, j))) / dx);

    coeff = _grid->coef(_x[i], _y[j]) / 2;
    if (_grid->get_xiy2_2(i, j) != 0)
      coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

    coeff = _grid->coef(_x[i], _y[j] - dy / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1.)) / dx);
    m = m + 1;
  }

  for (int j = 1; j < Ny - 1; ++j)
  {
    mat1[m] = _grid->ind(0, j) - 1;
    mat2[m] = _grid->ind(0, j) - 1;
    mat3[m] = 0;
    i = 0;

    coeff = _grid->coef(_x[0] + dx / 2, _y[j]);
    if (_grid->get_xix2_2(0, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(0, j) / 2 + _x[0] / 2, _y[j]);
    mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(0, j)) + _grid->get_xix2_2(0, j) * abs(_grid->get_xinter(0, j) - _x[0])) / dx);

    coeff = _grid->coef(_x[0] - dx / 2, _y[j]);
    mat3[m] = mat3[m] + coeff * (-1. / (dx * (1.)) / dx);

    coeff = _grid->coef(_x[0], _y[j] + dy / 2.);
    if (_grid->get_xiy2_2(0, j) != 0)
      coeff = _grid->coef(_x[0], _grid->get_yinter(0, j) / 2 + _y[j] / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, j)) + _grid->get_xiy2_2(0, j) * abs(_grid->get_yinter(0, j) - _y[j])) / dx);

    coeff = _grid->coef(_x[0], _y[j] - dy / 2);
    if (_grid->get_xiy2_2(0, j - 1) != 0)
      coeff = _grid->coef(_x[0], _grid->get_yinter(0, j - 1) / 2. + _y[j] / 2.);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, j - 1)) + _grid->get_xiy2_2(0, j - 1) * abs(_y[j] - _grid->get_yinter(0, j - 1))) / dx);
    m = m + 1;
  }

  for (int j = 1; j < Ny - 1; ++j)
  {
    int i = Nx - 1;
    mat1[m] = _grid->ind(i, j) - 1;
    mat2[m] = _grid->ind(i, j) - 1;
    mat3[m] = 0;

    coeff = _grid->coef(_x[i] + dx / 2, _y[j]);
    if (_grid->get_xix2_2(i, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx);

    coeff = _grid->coef(_x[i] - dx / 2, _y[j]);
    mat3[m] = mat3[m] + coeff * (-1. / (dx * (1.)) / dx);

    coeff = _grid->coef(_x[i], _y[j] + dy / 2.);
    if (_grid->get_xiy2_2(i, j) != 0)
      coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

    coeff = _grid->coef(_x[i], _y[j] - dy / 2);
    if (_grid->get_xiy2_2(i, j - 1) != 0)
      coeff = _grid->coef(_x[i], _grid->get_yinter(i, j - 1) / 2. + _y[j] / 2.);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j - 1)) + _grid->get_xiy2_2(i, j - 1) * abs(_y[j] - _grid->get_yinter(i, j - 1))) / dx);
    m = m + 1;
  }

  ///////////////////////////////////////////////////////////////////////////////

  for (int j = 1; j < Ny - 1; ++j)
  {
    for (int i = 1; i < Nx - 1; ++i)
    {

      mat1[m] = _grid->ind(i, j) - 1;
      mat2[m] = _grid->ind(i, j) - 1;
      mat3[m] = 0;

      coeff = (sigma(_grid->ind(i, j) - 1) + sigma(_grid->ind(i + 1, j) - 1)) / 2;
      // coeff = _grid->coef(_x[i] + dx / 2, _y[j]);
      if (_grid->get_xix2_2(i, j) != 0)
        coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
      mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx);

      coeff = (sigma(_grid->ind(i - 1, j) - 1) + sigma(_grid->ind(i, j) - 1)) / 2;
      // coeff = _grid->coef(_x[i] - dx / 2, _y[j]);
      if (_grid->get_xix2_2(i - 1, j) != 0)
        coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2 + _x[i] / 2, _y[j]);
      mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i - 1, j)) + _grid->get_xix2_2(i - 1, j) * abs(_x[i] - _grid->get_xinter(i - 1, j))) / dx);

      coeff = (sigma(_grid->ind(i, j) - 1) + sigma(_grid->ind(i, j + 1) - 1)) / 2;
      // coeff = _grid->coef(_x[i], _y[j] + dy / 2.);
      if (_grid->get_xiy2_2(i, j) != 0)
        coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
      mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

      coeff = (sigma(_grid->ind(i, j - 1) - 1) + sigma(_grid->ind(i, j) - 1)) / 2;
      // coeff = _grid->coef(_x[i], _y[j] - dy / 2);
      if (_grid->get_xiy2_2(i, j - 1) != 0)
        coeff = _grid->coef(_x[i], _grid->get_yinter(i, j - 1) / 2. + _y[j] / 2.);
      mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j - 1)) + _grid->get_xiy2_2(i, j - 1) * abs(_y[j] - _grid->get_yinter(i, j - 1))) / dx);

      m = m + 1;
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  for (int j = 0; j < Ny; ++j)
  {
    for (int i = 0; i < Nx; ++i)
    {

      if (i - 1 >= 0)
      {
        coeff = (sigma(_grid->ind(i - 1, j) - 1) + sigma(_grid->ind(i, j) - 1)) / 2;
        // coeff = _grid->coef(_x[i] - dx / 2., _y[j]);
        if (_grid->get_xix2_2(i - 1, j) != 0)
          coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2. + _x[i] / 2., _y[j]);
        mat1[m] = _grid->ind(i, j) - 1;
        mat2[m] = (1 - _grid->get_xix2_2(i - 1, j)) * _grid->ind(i - 1, j) + _grid->get_xix2_2(i - 1, j) * _grid->get_tx2(i - 1, j) - 1;
        mat3[m] = coeff * ((1 - _grid->get_xix2_2(i - 1, j)) / (dx * dx) + _grid->get_xix2_2(i - 1, j) / (abs(_x[i] - _grid->get_xinter(i - 1, j)) * dx));
        m = m + 1;
      }
      // i+1,j
      if (i + 1 <= Nx - 1)
      {
        coeff = (sigma(_grid->ind(i, j) - 1) + sigma(_grid->ind(i + 1, j) - 1)) / 2;
        // coeff = _grid->coef(_x[i] + dx / 2., _y[j]);
        if (_grid->get_xix2_2(i, j) != 0)
          coeff = _grid->coef(_grid->get_xinter(i, j) / 2. + _x[i] / 2, _y[j]);
        mat1[m] = _grid->ind(i, j) - 1;
        mat2[m] = ((1 - _grid->get_xix2_2(i, j)) * _grid->ind(i + 1, j) + _grid->get_xix2_2(i, j) * _grid->get_tx2(i, j)) - 1;
        mat3[m] = coeff * ((1 - _grid->get_xix2_2(i, j)) / (dx * dx) + _grid->get_xix2_2(i, j) / (abs(_grid->get_xinter(i, j) - _x[i]) * dx));
        m = m + 1;
      }
      // i,j-1
      if (j - 1 >= 0)
      {
        coeff = (sigma(_grid->ind(i, j - 1) - 1) + sigma(_grid->ind(i, j) - 1)) / 2;
        // coeff = _grid->coef(_x[i], _y[j] - dy / 2.);
        if (_grid->get_xiy2_2(i, j - 1) != 0)
          coeff = _grid->coef(_x[i], _grid->get_yinter(i, j - 1) / 2. + _y[j] / 2.);
        mat1[m] = _grid->ind(i, j) - 1;
        mat2[m] = (1 - _grid->get_xiy2_2(i, j - 1)) * _grid->ind(i, j - 1) + _grid->get_xiy2_2(i, j - 1) * _grid->get_ty2(i, j - 1) - 1;
        mat3[m] = coeff * ((1 - _grid->get_xiy2_2(i, j - 1)) / (dy * dy) + _grid->get_xiy2_2(i, j - 1) / (abs(_y[j] - _grid->get_yinter(i, j - 1)) * dy));
        m = m + 1;
      }
      // i,j+1
      if (j + 1 <= Nx - 1)
      {
        coeff = (sigma(_grid->ind(i, j) - 1) + sigma(_grid->ind(i, j + 1) - 1)) / 2;
        // coeff = _grid->coef(_x[i], _y[j] + dy / 2.);
        if (_grid->get_xiy2_2(i, j) != 0)
          coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2. + _y[j] / 2.);
        mat1[m] = _grid->ind(i, j) - 1;
        mat2[m] = ((1 - _grid->get_xiy2_2(i, j)) * _grid->ind(i, j + 1) + _grid->get_xiy2_2(i, j) * _grid->get_ty2(i, j)) - 1;
        mat3[m] = coeff * ((1 - _grid->get_xiy2_2(i, j)) / (dy * dy) + _grid->get_xiy2_2(i, j) / (abs(_grid->get_yinter(i, j) - _y[j]) * dy));
        m = m + 1;
      }
    }
  }

  for (int j = 1; j < Ny - 1; ++j)
  {
    for (int i = 1; i < Nx - 1; ++i)
    {

      if ((_grid->get_tx2(i, j)) != 0)
      {
        if (_grid->get_zone(_grid->get_tx2(i, j)) == 0)
          fluxx(i, j);
        else
          electrodefluxx(i, j, Z);
      }

      if ((_grid->get_ty2(i, j)) != 0)
      {
        if (_grid->get_zone(_grid->get_ty2(i, j)) == 0)
          fluxy(i, j);
        else
          electrodefluxy(i, j, Z);
      }
    }
  }

  Um.setConstant(0);
  for (int i = 1; i < Nx - 1; i++)
  {
    for (int j = 1; j < Ny - 1; j++)
    {

      if (_grid->get_tx2(i, j) != 0)
      {
        if (_grid->get_zone(_grid->get_tx2(i, j)) != 0)
        {
          Um(_grid->get_zone(_grid->get_tx2(i, j)) - 1) = Um(_grid->get_zone(_grid->get_tx2(i, j)) - 1) + dx * dy * _grid->Delta(i, j, alpha);
          mat1(m) = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i, j)) - 1;
          mat2(m) = _grid->ind(i, j) - 1;
          mat3[m] = -dx * dy * _grid->Delta(i, j, alpha);
          m = m + 1;
        }
      }
      else
      {
        if (_grid->get_tx2(i - 1, j) != 0)
        {
          if (_grid->get_zone(_grid->get_tx2(i - 1, j)) != 0)
          {
            Um(_grid->get_zone(_grid->get_tx2(i - 1, j)) - 1) = Um(_grid->get_zone(_grid->get_tx2(i - 1, j)) - 1) + dx * dy * _grid->Delta(i, j, alpha);
            mat1(m) = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i - 1, j)) - 1;
            mat2(m) = _grid->ind(i, j) - 1;
            mat3[m] = -dx * dy * _grid->Delta(i, j, alpha);
            m = m + 1;
          }
        }
        else
        {
          if (_grid->get_ty2(i, j) != 0)
          {
            if (_grid->get_zone(_grid->get_ty2(i, j)) != 0)
            {
              Um(_grid->get_zone(_grid->get_ty2(i, j)) - 1) = Um(_grid->get_zone(_grid->get_ty2(i, j)) - 1) + dx * dy * _grid->Delta(i, j, alpha);
              mat1(m) = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i, j)) - 1;
              mat2(m) = _grid->ind(i, j) - 1;
              mat3[m] = -dx * dy * _grid->Delta(i, j, alpha);
              m = m + 1;
            }
          }
          else
          {
            if (_grid->get_ty2(i, j - 1) != 0)
            {
              if (_grid->get_zone(_grid->get_ty2(i, j - 1)) != 0)
              {
                Um(_grid->get_zone(_grid->get_ty2(i, j - 1)) - 1) = Um(_grid->get_zone(_grid->get_ty2(i, j - 1)) - 1) + dx * dy * _grid->Delta(i, j, alpha);
                mat1(m) = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i, j - 1)) - 1;
                mat2(m) = _grid->ind(i, j) - 1;
                mat3[m] = -dx * dy * _grid->Delta(i, j, alpha);
                m = m + 1;
              }
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < Ne; i++)
  {
    mat1(m) = _grid->Get_nbsigne() + i;
    mat2(m) = _grid->Get_nbsigne() + i;
    mat3[m] = Um(i);
    if (i == 0)
      mat3[m] = mat3[m] + 0.0000000001;

    m = m + 1;
  }

  m = m - 1;

  matt1.resize(m + 1);
  matt2.resize(m + 1);
  matt3.resize(m + 1);

  for (int i = 0; i < m + 1; i++)
  {
    matt1[i] = mat1[i];
    matt2[i] = mat2[i];
    matt3[i] = mat3[i];
    trp.push_back(Trip(matt1[i], matt2[i], matt3[i]));
  }

  _lap_mat.setFromTriplets(trp.begin(), trp.end());
}

void EITAssembler::CLEAR()
{
  trp.clear();
  matt1.resize(0);
  matt2.resize(0);
  matt3.resize(0);
}

void EITAssembler::BuildSourceTerm(VectorXd &alpha, VectorXd &Imm, VectorXd &Z)
{
  _source_term.resize(Nx * Ny + _grid->Get_nbsigne() + Ne);
  _source_term.setConstant(0);

  Imm.resize(Ne);

  for (int i = 0; i < Ne; ++i)
  {
    _source_term[_grid->Get_nbsigne() + i + 1 - 1] = (1. / Z(i)) * Imm[i];
  }
}

void EITAssembler::BuildSourceTermMAN(VectorXd &alpha, VectorXd &Z)
{
  _source_term_man.resize(Nx * Ny + _grid->Get_nbsigne() + Ne);

  for (int j = 0; j < Ny; ++j)
  {
    for (int i = 0; i < Nx; ++i)
    {
      _source_term_man[_grid->ind(i, j) - 1] = _grid->source(_x[i], _y[j]);
    }
  }

  for (int j = 0; j < Ny; ++j)
  {
    for (int i = 0; i < Nx; ++i)
    {
      if (_grid->get_tx2(i, j) != 0)
      {
        if (_grid->get_zone(_grid->get_tx2(i, j)) == 0)
        {
          fix = (_grid->get_u(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_u(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
          fiy = (_grid->get_v(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_v(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
          fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
          fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
          fix = fiix;
          fiy = fiiy;
          _source_term_man[_grid->get_tx2(i, j) - 1] = _grid->val_du(_grid->get_xinter(i, j), _y[j], fix, fiy);
        }
        else
        {
          double elec_num = _grid->get_zone(_grid->get_tx2(i, j)) - 1;
          fix = (_grid->get_u(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_u(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
          fiy = (_grid->get_v(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_v(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
          fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
          fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
          fix = fiix;
          fiy = fiiy;
          _source_term_man[_grid->get_tx2(i, j) - 1] = Z[elec_num] * _grid->cond_inter(_grid->get_xinter(i, j), _y[j]) + _grid->val_du(_grid->get_xinter(i, j), _y[j], fix, fiy);
        }
      }
      if (_grid->get_ty2(i, j) != 0)
      {
        if (_grid->get_zone(_grid->get_ty2(i, j)) == 0)
        {
          fix = (_grid->get_u(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_u(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
          fiy = (_grid->get_v(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_v(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
          fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
          fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
          fix = fiix;
          fiy = fiiy;
          _source_term_man[_grid->get_ty2(i, j) - 1] = _grid->val_du(_x[i], _grid->get_yinter(i, j), fix, fiy); // ajouter les coef dee sigma
        }
        else
        {
          double elec_num = _grid->get_zone(_grid->get_ty2(i, j)) - 1;
          fix = (_grid->get_u(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_u(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
          fiy = (_grid->get_v(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_v(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
          fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
          fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
          fix = fiix;
          fiy = fiiy;
          _source_term_man[_grid->get_ty2(i, j) - 1] = Z[elec_num] * _grid->cond_inter(_x[i], _grid->get_yinter(i, j)) + _grid->val_du(_x[i], _grid->get_yinter(i, j), fix, fiy);
        }
      }
    }
  }

  Ummm.resize(Ne);
  Ummm.setConstant(0);

  for (int i = 1; i < Nx - 1; i++)
  {
    for (int j = 1; j < Ny - 1; j++)
    {

      if (_grid->get_tx2(i, j) != 0)
      {
        if (_grid->get_zone(_grid->get_tx2(i, j)) != 0)
        {
          Ummm(_grid->get_zone(_grid->get_tx2(i, j)) - 1) = Ummm(_grid->get_zone(_grid->get_tx2(i, j)) - 1) + _grid->cond_inter(_x[i], _y[j]) * dx * dy * _grid->Delta(i, j, alpha);
        }
      }
      else
      {
        if (_grid->get_tx2(i - 1, j) != 0)
        {
          if (_grid->get_zone(_grid->get_tx2(i - 1, j)) != 0)
          {
            Ummm(_grid->get_zone(_grid->get_tx2(i - 1, j)) - 1) = Ummm(_grid->get_zone(_grid->get_tx2(i - 1, j)) - 1) + _grid->cond_inter(_x[i], _y[j]) * dx * dy * _grid->Delta(i, j, alpha);
          }
        }
        else
        {
          if (_grid->get_ty2(i, j) != 0)
          {
            if (_grid->get_zone(_grid->get_ty2(i, j)) != 0)
            {
              Ummm(_grid->get_zone(_grid->get_ty2(i, j)) - 1) = Ummm(_grid->get_zone(_grid->get_ty2(i, j)) - 1) + _grid->cond_inter(_x[i], _y[j]) * dx * dy * _grid->Delta(i, j, alpha);
            }
          }
          else
          {
            if (_grid->get_ty2(i, j - 1) != 0)
            {
              if (_grid->get_zone(_grid->get_ty2(i, j - 1)) != 0)
              {
                Ummm(_grid->get_zone(_grid->get_ty2(i, j - 1)) - 1) = Ummm(_grid->get_zone(_grid->get_ty2(i, j - 1)) - 1) + _grid->cond_inter(_x[i], _y[j]) * dx * dy * _grid->Delta(i, j, alpha);
              }
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < Ne; ++i)
  {
    _source_term_man[_grid->Get_nbsigne() + i + 1 - 1] = 0 * Um(i) - (1. / 1) * Ummm(i);
  }
}

double EITAssembler::Compute_integral(MatrixXd &u,
                                      MatrixXd &w,
                                      VectorXd &alpha,
                                      VectorXd &angle_debut,
                                      VectorXd &angle_fin,
                                      int i,
                                      int l)
{

  VectorXd integ(Ne);

  integ.setConstant(0);

  for (int i = 1; i < Nx - 1; i++)
  {
    for (int j = 1; j < Ny - 1; j++)
    {

      if (_grid->get_tx2(i, j) != 0)
      {
        if (_grid->get_zone(_grid->get_tx2(i, j)) != 0)
        {
          double valu = u(l, _grid->ind(i, j)) - u(l, _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i, j)) - 1);
          double valw = w(l, _grid->ind(i, j)) - w(l, _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i, j)) - 1);
          integ(_grid->get_zone(_grid->get_tx2(i, j)) - 1) += valu * valw * dx * dy * _grid->Delta(i, j, alpha);
        }
      }
      else
      {
        if (_grid->get_tx2(i - 1, j) != 0)
        {
          if (_grid->get_zone(_grid->get_tx2(i - 1, j)) != 0)
          {
            double valu = u(l, _grid->ind(i, j)) - u(l, _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i - 1, j)) - 1);
            double valw = w(l, _grid->ind(i, j)) - w(l, _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i - 1, j)) - 1);
            integ(_grid->get_zone(_grid->get_tx2(i - 1, j)) - 1) += valu * valw * dx * dy * _grid->Delta(i, j, alpha);
          }
          else
          {
            if (_grid->get_ty2(i, j) != 0)
            {
              if (_grid->get_zone(_grid->get_ty2(i, j)) != 0)
              {
                double valu = u(l, _grid->ind(i, j)) - u(l, _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i, j)) - 1);
                double valw = w(l, _grid->ind(i, j)) - w(l, _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i, j)) - 1);
                integ(_grid->get_zone(_grid->get_ty2(i, j)) - 1) += valu * valw * dx * dy * _grid->Delta(i, j, alpha);
              }
            }
            else
            {
              if (_grid->get_ty2(i, j - 1) != 0)
              {
                if (_grid->get_zone(_grid->get_ty2(i, j - 1)) != 0)
                {
                  double valu = u(l, _grid->ind(i, j)) - u(l, _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i, j - 1)) - 1);
                  double valw = w(l, _grid->ind(i, j)) - w(l, _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i, j - 1)) - 1);
                  integ(_grid->get_zone(_grid->get_ty2(i, j - 1)) - 1) += valu * valw * dx * dy * _grid->Delta(i, j, alpha);
                }
              }
            }
          }
        }
      }
    }
  }

  return integ(i);
}

double EITAssembler::grad_normal_x(Eigen::MatrixXd &u, int i, int j, int l)
{

  pi = acos(-1.);

  fix = 0;
  fiy = 0;

  // calcul normale sur le point d'interface
  fix = (_grid->get_u(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_u(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
  fiy = (_grid->get_v(i, j) * (_x[i + 1] - _grid->get_xinter(i, j)) - _grid->get_v(i + 1, j) * (_x[i] - _grid->get_xinter(i, j))) / dx;
  fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
  fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
  fix = fiix;
  fiy = fiiy;

  // cout<<M << "\n"<<endl;

  Eigen::MatrixXd M, M_ind;
  M.resize(6, 2);
  M_ind.resize(6, 2);

  M(0, 0) = _x[i];
  M(0, 1) = _y[j + 1];
  M(1, 0) = _x[i + 1];
  M(1, 1) = _y[j + 1];
  M(2, 0) = _x[i];
  M(2, 1) = _y[j];
  M(3, 0) = _x[i + 1];
  M(3, 1) = _y[j];
  M(4, 0) = _x[i];
  M(4, 1) = _y[j - 1];
  M(5, 0) = _x[i + 1];
  M(5, 1) = _y[j - 1];

  M_ind(0, 0) = i;
  M_ind(0, 1) = j + 1;
  M_ind(1, 0) = i + 1;
  M_ind(1, 1) = j + 1;
  M_ind(2, 0) = i;
  M_ind(2, 1) = j;
  M_ind(3, 0) = i + 1;
  M_ind(3, 1) = j;
  M_ind(4, 0) = i;
  M_ind(4, 1) = j - 1;
  M_ind(5, 0) = i + 1;
  M_ind(5, 1) = j - 1;

  xj = _grid->get_xinter(i, j);
  yj = _y[j];
  _grid->discret_scheme(i, j, fix, fiy, xj, yj, M, M_ind);

  xk = _grid->get_xk(i, j);
  yk = _grid->get_yk(i, j);
  indk = _grid->get_indk(i, j);

  xi = _grid->get_xi(i, j);
  yi = _grid->get_yi(i, j);
  indi = _grid->get_indi(i, j);

  // calcul des coef du gradient (interpolation lineaire sur les 3 points)
  denom = (xj - xk) * (yj - yi) - (xj - xi) * (yj - yk);
  alphaj = (yk - yi) / denom;
  betaj = (xi - xk) / denom;

  denom = (xi - xk) * (yi - yj) - (xi - xj) * (yi - yk);
  alphai = (yk - yj) / denom;
  betai = (xj - xk) / denom;

  denom = (xk - xj) * (yk - yi) - (xk - xi) * (yk - yj);
  alphak = (yj - yi) / denom;
  betak = (xi - xj) / denom;

  return u(l, _grid->get_tx2(i, j) - 1) * (alphaj * fix + betaj * fiy) + u(l, indi - 1) * (alphai * fix + betai * fiy) + u(l, indk - 1) * (alphak * fix + betak * fiy);
}

double EITAssembler::grad_normal_y(Eigen::MatrixXd &u, int i, int j, int l)
{

  pi = acos(-1.);

  fix = 0;
  fiy = 0;

  // calcul normale sur le point d'interface
  fix = (_grid->get_u(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_u(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
  fiy = (_grid->get_v(i, j) * (_y[j + 1] - _grid->get_yinter(i, j)) - _grid->get_v(i, j + 1) * (_y[j] - _grid->get_yinter(i, j))) / dy;
  fiix = fix / sqrt(pow(fix, 2) + pow(fiy, 2));
  fiiy = fiy / sqrt(pow(fix, 2) + pow(fiy, 2));
  fix = fiix;
  fiy = fiiy;

  Eigen::MatrixXd M, M_ind;
  M.resize(6, 2);
  M_ind.resize(6, 2);

  M(0, 0) = _x[i - 1];
  M(0, 1) = _y[j + 1];
  M(1, 0) = _x[i - 1];
  M(1, 1) = _y[j];
  M(2, 0) = _x[i];
  M(2, 1) = _y[j + 1];
  M(3, 0) = _x[i];
  M(3, 1) = _y[j];
  M(4, 0) = _x[i + 1];
  M(4, 1) = _y[j + 1];
  M(5, 0) = _x[i + 1];
  M(5, 1) = _y[j];

  M_ind(0, 0) = i - 1;
  M_ind(0, 1) = j + 1;
  M_ind(1, 0) = i - 1;
  M_ind(1, 1) = j;
  M_ind(2, 0) = i;
  M_ind(2, 1) = j + 1;
  M_ind(3, 0) = i;
  M_ind(3, 1) = j;
  M_ind(4, 0) = i + 1;
  M_ind(4, 1) = j + 1;
  M_ind(5, 0) = i + 1;
  M_ind(5, 1) = j;

  xj = _x[i];
  yj = _grid->get_yinter(i, j);
  _grid->discret_scheme(i, j, fix, fiy, xj, yj, M, M_ind);

  xk = _grid->get_xk(i, j);
  yk = _grid->get_yk(i, j);
  indk = _grid->get_indk(i, j);

  xi = _grid->get_xi(i, j);
  yi = _grid->get_yi(i, j);
  indi = _grid->get_indi(i, j);

  // calcul des coef du gradient (interpolation lineaire sur les 3 points)
  denom = (xj - xk) * (yj - yi) - (xj - xi) * (yj - yk);
  alphaj = (yk - yi) / denom;
  betaj = (xi - xk) / denom;

  denom = (xi - xk) * (yi - yj) - (xi - xj) * (yi - yk);
  alphai = (yk - yj) / denom;
  betai = (xj - xk) / denom;

  denom = (xk - xj) * (yk - yi) - (xk - xi) * (yk - yj);
  alphak = (yj - yi) / denom;
  betak = (xi - xj) / denom;

  return u(l, _grid->get_ty2(i, j) - 1) * (alphaj * fix + betaj * fiy) + u(l, indi - 1) * (alphai * fix + betai * fiy) + u(l, indk - 1) * (alphak * fix + betak * fiy);
}

void EITAssembler::Buildsourcedir1(uint l, MatrixXd u, MatrixXd w, VectorXd &alpha)
{

  // _grid->BuildT(alpha);
  // _grid->inter(alpha);

  // u.resize(Ne, Nx * Ny + _grid->Get_nbsigne() + Ne);
  // w.resize(Ne, Nx * Ny + _grid->Get_nbsigne() + Ne);
  dsigma.resize(Nx, Ny);
  double uup, uum, wp, wm;
  double xp, xm, yp, ym;

  for (int i = 0; i < Nx; ++i)
  {
    for (int j = 0; j < Ny; ++j)
    {
      dsigma(i, j) = 0;
    }
  }

  for (int i = 0; i < Nx; ++i)
  {
    for (int j = 0; j < Ny; ++j)
    {

      if (_grid->Levelset(i, j, alpha) <= 0.)
      {

        uup = u(l, _grid->ind(i + 1, j) - 1);
        wp = w(l, _grid->ind(i + 1, j) - 1);
        xp = _x[i + 1];
        if (_grid->get_tx2(i, j) != 0)
        {
          uup = u(l, _grid->get_tx2(i, j) - 1);
          wp = w(l, _grid->get_tx2(i, j) - 1);
          xp = _grid->get_xinter(i, j);
        }

        uum = u(l, _grid->ind(i - 1, j) - 1);
        wm = w(l, _grid->ind(i - 1, j) - 1);
        xm = _x[i - 1];
        if (_grid->get_tx2(i - 1, j) != 0)
        {
          uum = u(l, _grid->get_tx2(i - 1, j) - 1);
          wm = w(l, _grid->get_tx2(i - 1, j) - 1);
          xm = _grid->get_xinter(i - 1, j);
        }

        dsigma(i, j) = dsigma(i, j) + ((uup - uum) / (xp - xm) * (wp - wm) / (xp - xm));

        uup = u(l, _grid->ind(i, j + 1) - 1);
        wp = w(l, _grid->ind(i, j + 1) - 1);
        yp = _y[j + 1];
        if (_grid->get_ty2(i, j) != 0)
        {
          uup = u(l, _grid->get_ty2(i, j) - 1);
          wp = w(l, _grid->get_ty2(i, j) - 1);
          yp = _grid->get_yinter(i, j);
        }

        uum = u(l, _grid->ind(i, j - 1) - 1);
        wm = w(l, _grid->ind(i, j - 1) - 1);
        ym = _y[j - 1];
        if (_grid->get_ty2(i, j - 1) != 0)
        {
          uum = u(l, _grid->get_ty2(i, j - 1) - 1);
          wm = w(l, _grid->get_ty2(i, j - 1) - 1);
          ym = _grid->get_yinter(i, j - 1);
        }

        dsigma(i, j) = dsigma(i, j) + ((uup - uum) / (yp - ym)) * ((wp - wm) / (yp - ym));
      }
    }
  }
}

void EITAssembler::Buildsourcedir2(Eigen::VectorXd u, Eigen::VectorXd v, Eigen::VectorXd alpha)
{

  u.resize(Nx * Ny + _grid->Get_nbsigne() + Ne);
  v.resize(Nx * Ny + _grid->Get_nbsigne() + Ne);
  lap_dir.resize(Nx, Ny);

  for (int i = 0; i < Nx; ++i)
  {
    for (int j = 0; j < Ny; ++j)
    {
      lap_dir(i, j) = 0;
    }
  }

  for (int i = 0; i < Nx; ++i)
  {
    for (int j = 0; j < Ny; ++j)
    {
      if (_grid->Levelset(i, j, alpha) <= 0.)
      {
        lap_dir(i, j) = (((u(_grid->ind(i + 1, j) - 1) - v(_grid->ind(i + 1, j) - 1)) -
                          2 * (u(_grid->ind(i, j) - 1) - v(_grid->ind(i, j) - 1)) +
                          (u(_grid->ind(i - 1, j) - 1) - v(_grid->ind(i - 1, j) - 1))) /
                             pow(dx, 2) +
                         ((u(_grid->ind(i, j + 1) - 1) - v(_grid->ind(i, j + 1) - 1)) - 2 * (u(_grid->ind(i, j) - 1) - v(_grid->ind(i, j) - 1)) +
                          (u(_grid->ind(i, j - 1) - 1) - v(_grid->ind(i, j - 1) - 1))) /
                             pow(dy, 2));
      }
    }
  }
}



VectorXd EITAssembler::grad_of_vec(VectorXd alpha, VectorXd vec)
{

  VectorXd grad;

  grad = vec;

  grad.setConstant(0);
  double uup, uum, xp, xm, yp, ym;

  for (int i = 0; i < _grid->getNbPoints(0); i++)
  {
    for (int j = 0; j < _grid->getNbPoints(1); j++)
    {
      if (_grid->Levelset(i, j, alpha) <= 0)
      {

        double uup = vec(_grid->ind(i + 1, j) - 1);
        double xp = _x[i + 1];
        if (_grid->get_tx2(i, j) != 0)
        {
          uup = vec(_grid->get_tx2(i, j) - 1);
          xp = _grid->get_xinter(i, j);
        }

        double uum = vec(_grid->ind(i - 1, j) - 1);
        double xm = _x[i - 1];
        if (_grid->get_tx2(i - 1, j) != 0)
        {
          uum = vec(_grid->get_tx2(i - 1, j) - 1);
          xm = _grid->get_xinter(i - 1, j);
        }
        grad(_grid->ind(i, j)-1) += ((uup - uum) / (xp - xm));

        uup = vec(_grid->ind(i, j + 1) - 1);
        yp = _y[j + 1];
        if (_grid->get_ty2(i, j) != 0)
        {
          uup = vec(_grid->get_ty2(i, j) - 1);
          yp = _grid->get_yinter(i, j);
        }

        uum = vec(_grid->ind(i, j - 1) - 1);
        ym = _y[j - 1];
        if (_grid->get_ty2(i, j - 1) != 0)
        {
          uum = vec(_grid->get_ty2(i, j - 1) - 1);
          ym = _grid->get_yinter(i, j - 1);
        }
        grad(_grid->ind(i, j)-1) +=  ((uup - uum) / (yp - ym));
      }
    }
  }

  return grad;
}

VectorXd EITAssembler::Build_dirVT_G(VectorXd alpha, VectorXd vec)
{

  double eta = 0.001, den, nom;
  VectorXd G;

  nom = vec.cwiseAbs().dot(vec.cwiseAbs()) + 2 * eta;
  den = vec.cwiseAbs().dot(vec.cwiseAbs()) + eta;

  G = (0.001 * 0.5 * (nom) / (pow(den, 3 / 2))) * vec;

  return G;
}

void EITAssembler::BuildDirMatrix(VectorXd alpha, VectorXd angle_debut, VectorXd angle_fin, Eigen::VectorXd sigma)
{

  // _grid->BuildT(alpha);
  // _grid->inter(alpha);
  // _grid->DLevelset(alpha);
  // _grid->Zonelectrode(alpha, angle_debut, angle_fin);

  _lap_dir.resize(Nx * Ny + _grid->Get_nbsigne() + Ne, Nx * Ny + _grid->Get_nbsigne() + Ne);
  sigma.resize(Nx * Ny + _grid->Get_nbsigne() + Ne);

  double coeff = 0;
  m = 0;

  // Premier cas:: i=0; j=0.
  mat1[m] = _grid->ind(0, 0) - 1;
  mat2[m] = _grid->ind(0, 0) - 1;
  mat3[m] = 0;

  coeff = _grid->coef(_x[0] + dx / 2, _y[0]);
  if (_grid->get_xix2_2(0, 0) != 0)
    coeff = _grid->coef(_grid->get_xinter(0, 0) / 2 + _x[0] / 2, _y[0]);
  mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(0, 0)) + _grid->get_xix2_2(0, 0) * abs(_grid->get_xinter(0, 0) - _x[0])) / dx) - 1;

  coeff = _grid->coef(_x[0], _y[0]) / 2;
  mat3[m] = mat3[m] + coeff * (-1. / (dx * (1.)) / dx);

  coeff = _grid->coef(_x[0], _y[0] + dy / 2.);
  if (_grid->get_xiy2_2(0, 0) != 0)
    coeff = _grid->coef(_x[0], _grid->get_yinter(0, 0) / 2 + _y[0] / 2);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, 0)) + _grid->get_xiy2_2(0, 0) * abs(_grid->get_yinter(0, 0) - _y[0])) / dx);

  coeff = _grid->coef(_x[0], _y[0]) / 2;
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1.)) / dx);

  m = m + 1;

  // i=0, j=Ny-1
  mat1[m] = _grid->ind(0, Ny - 1) - 1;
  mat2[m] = _grid->ind(0, Ny - 1) - 1;
  mat3[m] = 0;

  coeff = _grid->coef(_x[0] + dx / 2, _y[Ny - 1]);
  if (_grid->get_xix2_2(0, Ny - 1) != 0)
    coeff = _grid->coef(_grid->get_xinter(0, Ny - 1) / 2 + _x[0] / 2, _y[Ny - 1]);
  mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(0, Ny - 1)) + _grid->get_xix2_2(0, Ny - 1) * abs(_grid->get_xinter(0, Ny - 1) - _x[0])) / dx) - 1;

  coeff = _grid->coef(_x[0] - dx / 2, _y[Ny - 1]);
  mat3[m] = mat3[m] + coeff * (-1. / (dx * (1.)) / dx);

  coeff = _grid->coef(_x[0], _y[Ny - 1]) / 2;
  if (_grid->get_xiy2_2(0, Ny - 1) != 0)
    coeff = _grid->coef(_x[0], _grid->get_yinter(0, Ny - 1) / 2 + _y[Ny - 1] / 2);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, Ny - 1)) + _grid->get_xiy2_2(0, Ny - 1) * abs(_grid->get_yinter(0, Ny - 1) - _y[Ny - 1])) / dx);

  coeff = _grid->coef(_x[0], _y[Ny - 1] - dy / 2);
  if (_grid->get_xiy2_2(0, Ny - 1 - 1) != 0)
    coeff = _grid->coef(_x[0], _grid->get_yinter(0, Ny - 1 - 1) / 2. + _y[Ny - 1] / 2.);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, Ny - 1 - 1)) + _grid->get_xiy2_2(0, Ny - 1 - 1) * abs(_y[Ny - 1] - _grid->get_yinter(0, Ny - 1 - 1))) / dx);

  m = m + 1;

  // i=Nx-1, j=0
  mat1[m] = _grid->ind(Nx - 1, 0) - 1;
  mat2[m] = _grid->ind(Nx - 1, 0) - 1;
  mat3[m] = 0;

  coeff = _grid->coef(_x[Nx - 1] + dx / 2, _y[0]);
  if (_grid->get_xix2_2(Nx - 1, 0) != 0)
    coeff = _grid->coef(_grid->get_xinter(Nx - 1, 0) / 2 + _x[Nx - 1] / 2, _y[0]);
  mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(Nx - 1, 0)) + _grid->get_xix2_2(Nx - 1, 0) * abs(_grid->get_xinter(Nx - 1, 0) - _x[Nx - 1])) / dx) - 1;

  coeff = _grid->coef(_x[Nx - 1] - dx / 2, _y[0]);
  if (_grid->get_xix2_2(Nx - 1 - 1, 0) != 0)
    coeff = _grid->coef(_grid->get_xinter(Nx - 1 - 1, 0) / 2 + _x[Nx - 1] / 2, _y[0]);
  mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(Nx - 1 - 1, 0)) + _grid->get_xix2_2(Nx - 1 - 1, 0) * abs(_x[Nx - 1] - _grid->get_xinter(Nx - 1 - 1, 0))) / dx);

  coeff = _grid->coef(_x[Nx - 1], _y[0] + dy / 2.);
  if (_grid->get_xiy2_2(Nx - 1, 0) != 0)
    coeff = _grid->coef(_x[Nx - 1], _grid->get_yinter(Nx - 1, 0) / 2 + _y[0] / 2);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(Nx - 1, 0)) + _grid->get_xiy2_2(Nx - 1, 0) * abs(_grid->get_yinter(Nx - 1, 0) - _y[0])) / dx);

  coeff = _grid->coef(_x[Nx - 1], _y[0]) / 2;
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1.)) / dx);

  m = m + 1;

  // i=Nx-1, j=Ny-1
  mat1[m] = _grid->ind(Nx - 1, Ny - 1) - 1;
  mat2[m] = _grid->ind(Nx - 1, Ny - 1) - 1;
  mat3[m] = 0;

  int i = Nx - 1;
  int j = Ny - 1;

  coeff = _grid->coef(_x[i], _y[j]) / 2;
  if (_grid->get_xix2_2(i, j) != 0)
    coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
  mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx) - 1;

  coeff = _grid->coef(_x[i] - dx / 2, _y[j]);
  if (_grid->get_xix2_2(i - 1, j) != 0)
    coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2 + _x[i] / 2, _y[j]);
  mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i - 1, j)) + _grid->get_xix2_2(i - 1, j) * abs(_x[i] - _grid->get_xinter(i - 1, j))) / dx);

  coeff = _grid->coef(_x[i], _y[j]) / 2;
  if (_grid->get_xiy2_2(i, j) != 0)
    coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

  coeff = _grid->coef(_x[i], _y[j] - dy / 2);
  if (_grid->get_xiy2_2(i, j - 1) != 0)
    coeff = _grid->coef(_x[i], _grid->get_yinter(i, j - 1) / 2. + _y[j] / 2.);
  mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j - 1)) + _grid->get_xiy2_2(i, j - 1) * abs(_y[j] - _grid->get_yinter(i, j - 1))) / dx);

  m = m + 1;

  for (int i = 1; i < Nx - 1; ++i)
  {

    int j = 0;
    mat1[m] = _grid->ind(i, j) - 1;
    mat2[m] = _grid->ind(i, j) - 1;
    mat3[m] = 0;

    coeff = _grid->coef(_x[i] + dx / 2, _y[j]);
    if (_grid->get_xix2_2(i, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx) - 1;

    coeff = _grid->coef(_x[i] - dx / 2, _y[j]);
    if (_grid->get_xix2_2(i - 1, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i - 1, j)) + _grid->get_xix2_2(i - 1, j) * abs(_x[i] - _grid->get_xinter(i - 1, j))) / dx);

    coeff = _grid->coef(_x[i], _y[j] + dy / 2.);
    if (_grid->get_xiy2_2(i, j) != 0)
      coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

    coeff = _grid->coef(_x[i], _y[j]) / 2;
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1.)) / dx);

    m = m + 1;
  }

  for (int i = 1; i < Nx - 1; ++i)
  {

    int j = Ny - 1;
    mat1[m] = _grid->ind(i, j) - 1;
    mat2[m] = _grid->ind(i, j) - 1;
    mat3[m] = 0;

    coeff = _grid->coef(_x[i] + dx / 2, _y[j]);
    if (_grid->get_xix2_2(i, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx) - 1;

    coeff = _grid->coef(_x[i] - dx / 2, _y[j]);
    if (_grid->get_xix2_2(i - 1, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i - 1, j)) + _grid->get_xix2_2(i - 1, j) * abs(_x[i] - _grid->get_xinter(i - 1, j))) / dx);

    coeff = _grid->coef(_x[i], _y[j]) / 2;
    if (_grid->get_xiy2_2(i, j) != 0)
      coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

    coeff = _grid->coef(_x[i], _y[j] - dy / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1.)) / dx);

    m = m + 1;
  }

  for (int j = 1; j < Ny - 1; ++j)
  {

    mat1[m] = _grid->ind(0, j) - 1;
    mat2[m] = _grid->ind(0, j) - 1;
    mat3[m] = 0;

    coeff = _grid->coef(_x[0] + dx / 2, _y[j]);
    if (_grid->get_xix2_2(0, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(0, j) / 2 + _x[0] / 2, _y[j]);
    mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(0, j)) + _grid->get_xix2_2(0, j) * abs(_grid->get_xinter(0, j) - _x[0])) / dx) - 1;

    coeff = _grid->coef(_x[0] - dx / 2, _y[j]);
    mat3[m] = mat3[m] + coeff * (-1. / (dx * (1.)) / dx);

    coeff = _grid->coef(_x[0], _y[j] + dy / 2.);
    if (_grid->get_xiy2_2(0, j) != 0)
      coeff = _grid->coef(_x[0], _grid->get_yinter(0, j) / 2 + _y[j] / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, j)) + _grid->get_xiy2_2(0, j) * abs(_grid->get_yinter(0, j) - _y[j])) / dx);

    coeff = _grid->coef(_x[0], _y[j] - dy / 2);
    if (_grid->get_xiy2_2(0, j - 1) != 0)
      coeff = _grid->coef(_x[0], _grid->get_yinter(0, j - 1) / 2. + _y[j] / 2.);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(0, j - 1)) + _grid->get_xiy2_2(0, j - 1) * abs(_y[j] - _grid->get_yinter(0, j - 1))) / dx);

    m = m + 1;
  }

  for (int j = 1; j < Ny - 1; ++j)
  {

    int i = Nx - 1;
    mat1[m] = _grid->ind(i, j) - 1;
    mat2[m] = _grid->ind(i, j) - 1;
    mat3[m] = 0;

    coeff = _grid->coef(_x[i] + dx / 2, _y[j]);
    if (_grid->get_xix2_2(i, j) != 0)
      coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
    mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx) - 1;

    coeff = _grid->coef(_x[i] - dx / 2, _y[j]);
    mat3[m] = mat3[m] + coeff * (-1. / (dx * (1.)) / dx);

    coeff = _grid->coef(_x[i], _y[j] + dy / 2.);
    if (_grid->get_xiy2_2(i, j) != 0)
      coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

    coeff = _grid->coef(_x[i], _y[j] - dy / 2);
    if (_grid->get_xiy2_2(i, j - 1) != 0)
      coeff = _grid->coef(_x[i], _grid->get_yinter(i, j - 1) / 2. + _y[j] / 2.);
    mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j - 1)) + _grid->get_xiy2_2(i, j - 1) * abs(_y[j] - _grid->get_yinter(i, j - 1))) / dx);

    m = m + 1;
  }

  ///////////////////////////////////////////////////////////////////////////////
  for (int j = 1; j < Ny - 1; ++j)
  {
    for (int i = 1; i < Nx - 1; ++i)
    {

      mat1[m] = _grid->ind(i, j) - 1;
      mat2[m] = _grid->ind(i, j) - 1;
      mat3[m] = 0;

      coeff = (sigma(_grid->ind(i, j) - 1) + sigma(_grid->ind(i + 1, j) - 1)) / 2;

      if (_grid->get_xix2_2(i, j) != 0)
        coeff = _grid->coef(_grid->get_xinter(i, j) / 2 + _x[i] / 2, _y[j]);
      mat3[m] = coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i, j)) + _grid->get_xix2_2(i, j) * abs(_grid->get_xinter(i, j) - _x[i])) / dx) - 1;

      coeff = (sigma(_grid->ind(i - 1, j) - 1) + sigma(_grid->ind(i, j) - 1)) / 2;

      if (_grid->get_xix2_2(i - 1, j) != 0)
        coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2 + _x[i] / 2, _y[j]);
      mat3[m] = mat3[m] + coeff * (-1. / (dx * (1. - _grid->get_xix2_2(i - 1, j)) + _grid->get_xix2_2(i - 1, j) * abs(_x[i] - _grid->get_xinter(i - 1, j))) / dx);

      coeff = (sigma(_grid->ind(i, j) - 1) + sigma(_grid->ind(i, j + 1) - 1)) / 2;

      if (_grid->get_xiy2_2(i, j) != 0)
        coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2 + _y[j] / 2);
      mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j)) + _grid->get_xiy2_2(i, j) * abs(_grid->get_yinter(i, j) - _y[j])) / dx);

      coeff = (sigma(_grid->ind(i, j - 1) - 1) + sigma(_grid->ind(i, j) - 1)) / 2;

      if (_grid->get_xiy2_2(i, j - 1) != 0)
        coeff = _grid->coef(_x[i], _grid->get_yinter(i, j - 1) / 2. + _y[j] / 2.);
      mat3[m] = mat3[m] + coeff * (-1. / (dy * (1. - _grid->get_xiy2_2(i, j - 1)) + _grid->get_xiy2_2(i, j - 1) * abs(_y[j] - _grid->get_yinter(i, j - 1))) / dx);

      m = m + 1;
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  for (int j = 0; j < Ny; ++j)
  {
    for (int i = 0; i < Nx; ++i)
    {

      if (i - 1 >= 0)
      {
        coeff = (sigma(_grid->ind(i - 1, j) - 1) + sigma(_grid->ind(i, j) - 1)) / 2;
        if (_grid->get_xix2_2(i - 1, j) != 0)
          coeff = _grid->coef(_grid->get_xinter(i - 1, j) / 2. + _x[i] / 2., _y[j]);
        mat1[m] = _grid->ind(i, j) - 1;
        mat2[m] = (1 - _grid->get_xix2_2(i - 1, j)) * _grid->ind(i - 1, j) + _grid->get_xix2_2(i - 1, j) * _grid->get_tx2(i - 1, j) - 1;
        mat3[m] = coeff * ((1 - _grid->get_xix2_2(i - 1, j)) / (dx * dx) + _grid->get_xix2_2(i - 1, j) / (abs(_x[i] - _grid->get_xinter(i - 1, j)) * dx));
        m = m + 1;
      }
      // i+1,j
      if (i + 1 <= Nx - 1)
      {
        coeff = (sigma(_grid->ind(i, j) - 1) + sigma(_grid->ind(i + 1, j) - 1)) / 2;
        if (_grid->get_xix2_2(i, j) != 0)
          coeff = _grid->coef(_grid->get_xinter(i, j) / 2. + _x[i] / 2, _y[j]);
        mat1[m] = _grid->ind(i, j) - 1;
        mat2[m] = ((1 - _grid->get_xix2_2(i, j)) * _grid->ind(i + 1, j) + _grid->get_xix2_2(i, j) * _grid->get_tx2(i, j)) - 1;
        mat3[m] = coeff * ((1 - _grid->get_xix2_2(i, j)) / (dx * dx) + _grid->get_xix2_2(i, j) / (abs(_grid->get_xinter(i, j) - _x[i]) * dx));
        m = m + 1;
      }
      // i,j-1
      if (j - 1 >= 0)
      {
        coeff = (sigma(_grid->ind(i, j - 1) - 1) + sigma(_grid->ind(i, j) - 1)) / 2;
        if (_grid->get_xiy2_2(i, j - 1) != 0)
          coeff = _grid->coef(_x[i], _grid->get_yinter(i, j - 1) / 2. + _y[j] / 2.);
        mat1[m] = _grid->ind(i, j) - 1;
        mat2[m] = (1 - _grid->get_xiy2_2(i, j - 1)) * _grid->ind(i, j - 1) + _grid->get_xiy2_2(i, j - 1) * _grid->get_ty2(i, j - 1) - 1;
        mat3[m] = coeff * ((1 - _grid->get_xiy2_2(i, j - 1)) / (dy * dy) + _grid->get_xiy2_2(i, j - 1) / (abs(_y[j] - _grid->get_yinter(i, j - 1)) * dy));
        m = m + 1;
      }
      // i,j+1
      if (j + 1 <= Nx - 1)
      {
        coeff = (sigma(_grid->ind(i, j) - 1) + sigma(_grid->ind(i, j + 1) - 1)) / 2;
        if (_grid->get_xiy2_2(i, j) != 0)
          coeff = _grid->coef(_x[i], _grid->get_yinter(i, j) / 2. + _y[j] / 2.);
        mat1[m] = _grid->ind(i, j) - 1;
        mat2[m] = ((1 - _grid->get_xiy2_2(i, j)) * _grid->ind(i, j + 1) + _grid->get_xiy2_2(i, j) * _grid->get_ty2(i, j)) - 1;
        mat3[m] = coeff * ((1 - _grid->get_xiy2_2(i, j)) / (dy * dy) + _grid->get_xiy2_2(i, j) / (abs(_grid->get_yinter(i, j) - _y[j]) * dy));
        m = m + 1;
      }
    }
  }

  for (int j = 1; j < Ny - 1; ++j)
  {
    for (int i = 1; i < Nx - 1; ++i)
    {

      // s'il y a un point d'interface entre i et i+1
      if ((_grid->get_tx2(i, j)) != 0)
      {
        mat1[m] = _grid->get_tx2(i, j) - 1;
        mat2[m] = _grid->get_tx2(i, j) - 1;
        mat3[m] = 1;
        m = m + 1;
      }

      // s'il y a un point d'interface entre i et i+1
      if ((_grid->get_ty2(i, j)) != 0)
      {
        mat1[m] = _grid->get_ty2(i, j) - 1;
        mat2[m] = _grid->get_ty2(i, j) - 1;
        mat3[m] = 1;
        m = m + 1;
      }
    }
  }

  for (int i = 0; i < Ne; i++)
  {
    mat1(m) = _grid->Get_nbsigne() + i;
    mat2(m) = _grid->Get_nbsigne() + i;
    mat3[m] = 1;
    if (i == 0)
    {
      mat3[m] = mat3[m] + 0.0000000001;
    }
    m = m + 1;
  }

  m = m - 1;

  matt1.resize(m + 1);
  matt2.resize(m + 1);
  matt3.resize(m + 1);

  for (int i = 0; i < m + 1; i++)
  {
    matt1[i] = mat1[i];
    matt2[i] = mat2[i];
    matt3[i] = mat3[i];
    trp.push_back(Trip(matt1[i], matt2[i], matt3[i]));
  }

  _lap_dir.setFromTriplets(trp.begin(), trp.end());
}
