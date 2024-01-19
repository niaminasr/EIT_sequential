#include "grid.hpp"
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <cmath>
#include "Eigen/Eigen/Eigen"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

using namespace std;
using namespace Eigen;

Grid::Grid(Parameters &p)
{

  _N[0] = p.getInt("Nx");
  _N[1] = p.getInt("Ny");
  _N[2] = p.getInt("Ne");

  if (_N[0] < 3 or _N[1] < 3 or _N[2] < 1)
  {
    std::cerr << "Grid : Insufficient or negative number of points/ number of electrodes  : x " << _N[0] << " y " << _N[1] << "Ne" << _N[2] << std::endl;
    abort();
  }

  squaresize = 2.;

  _h[0] = (squaresize * 2) / (long double)(_N[0]);
  _h[1] = (squaresize * 2) / (long double)(_N[1]);

  _x.resize(_N[0]);
  for (unsigned short i = 0; i < ((unsigned short)(_N[0])); ++i)
    _x[i] = ((i + 1) - 0.5) * _h[0] - squaresize;

  _y.resize(_N[1]);
  for (unsigned short i = 0; i < ((unsigned short)(_N[1])); ++i)
    _y[i] = ((i + 1) - 0.5) * _h[1] - squaresize;
}

long double Grid::getSpacing(uint dim)
{

  if (dim != 0 and dim != 1)
  {
    std::cerr << "Grid::getSpacing: wrong dimension " << dim << std::endl;
    abort();
  }
  return _h[dim];
}

uint Grid::getNbPoints(uint dim)
{

  if (dim != 0 and dim != 1 and dim != 2)
  {
    std::cerr << "Grid::getNbPoints: wrong dimension " << dim << std::endl;
    abort();
  }
  return _N[dim];
}

void Grid::index2D(uint row, uint &i, uint &j)
{

  if (row >= _N[0] * _N[1])
  {
    std::cerr << "Grid::index2D: given row index (" << row << ")"
              << " is not in matrix of size " << _N[0] * _N[1] << std::endl;
    abort();
  }
  j = row / _N[0];
  i = row % _N[0];
}

void Grid::coords(uint i, uint j, double &x, double &y)
{
  x = i * _h[0];
  y = j * _h[1];
}

double Grid::rayon(double theta, VectorXd alpha)
{
  int n_alpha = alpha.size();
  int n = (n_alpha - 1) / 2;
  int i;
  double r_c, r_s;

  r_c = alpha(0);
  r_s = 0;
  for (int i = 1; i <= n; i++)
  {
    r_c = r_c + (alpha(i) * cos(i * theta));
    r_s = r_s + (alpha(n + i) * sin(i * theta));
  }

  return r_c + r_s;
}

long double Grid::rprime(double theta, VectorXd alpha)
{
  int n_alpha = alpha.size();
  int n = (n_alpha - 1) / 2;
  int i;
  double r_c, r_s;

  r_c = 0;
  r_s = 0;
  for (int i = 1; i <= n; i++)
  {
    r_c = r_c - (i * (alpha(i) * sin(i * theta)));
    r_s = r_s + (i * (alpha(n + i) * cos(i * theta)));
  }

  return r_c + r_s;
}

long double Grid::rdoubleprime(double theta, VectorXd alpha)
{
  int n_alpha = alpha.size();
  int n = (n_alpha - 1) / 2;
  int i;
  double r_c, r_s;

  r_c = 0;
  r_s = 0;
  for (int i = 1; i <= n; i++)
  {
    r_c = r_c - ((pow(i, 2)) * (alpha(i) * cos(i * theta)));
    r_s = r_s - ((pow(i, 2)) * (alpha(n + i) * sin(i * theta)));
  }
  return r_c + r_s;
}

double Grid::angle(double x, double y)
{
  return atan2(y, x);
}

double Grid::rho(double theta, VectorXd alpha)
{
  return sqrt(pow(rayon(theta, alpha), 2) + pow(rprime(theta, alpha), 2));
}

long double Grid::normalvec_x(double theta, VectorXd alpha)
{
  return (1 / rho(theta, alpha)) * (rayon(theta, alpha) * cos(theta) +
                                    rprime(theta, alpha) * sin(theta));
}

long double Grid::normalvec_y(double theta, VectorXd alpha)
{
  return (1 / rho(theta, alpha)) * (rayon(theta, alpha) * sin(theta) -
                                    rprime(theta, alpha) * cos(theta));
}

long double Grid::tangentvec_x(double theta, VectorXd alpha)
{
  return (1 / rho(theta, alpha)) * (rprime(theta, alpha) * cos(theta) -
                                    rayon(theta, alpha) * sin(theta));
}

long double Grid::tangentvec_y(double theta, VectorXd alpha)
{
  return (1 / rho(theta, alpha)) * (rprime(theta, alpha) * sin(theta) -
                                    rayon(theta, alpha) * cos(theta));
}

double Grid::funcI(double theta, VectorXd alpha)
{

  return 2 * exp(pow(rayon(theta, alpha), 2)) *
         ((pow(rayon(theta, alpha), 2)) / rho(theta, alpha)) * rho(theta, alpha);
}

double Grid::fintegral(double L, double theta_moin, double theta, Eigen::VectorXd alpha)
{
  int l, n;
  double d_the, Fint;

  n = 1000;
  d_the = (theta - theta_moin) / n;
  Fint = L;

  for (int i = 0; i < n; i++)
  {
    Fint = Fint - (d_the * rho(theta_moin + (i + 1 / 2) * d_the, alpha));
  }

  return Fint;
}

double Grid::angleplus(double L, double theta_moin, VectorXd alpha)
{
  double theta_0, theta1, epsilon;
  int l;
  double pi = acos(-1.);

  epsilon = 0.0001;
  theta_0 = theta_moin - L / 0.5;
  l = 0;

  while (abs(fintegral(L, theta_moin, theta_0, alpha)) >= epsilon)
  {
    theta1 = theta_0 + (fintegral(L, theta_moin, theta_0, alpha) / rho(theta_0, alpha));
    theta_0 = theta1;
    l = l + 1;
  }

  return theta1;
}

long double Grid::Levelset(uint i, uint j, Eigen::VectorXd alpha)
{
  double yc(0.), xc(0.), r;
  double r0(0.5);

  double xx, yy;

  xx = ((float(i) + 1) - 0.5) * _h[0] - squaresize;
  yy = ((float(j) + 1) - 0.5) * _h[1] - squaresize;

  r = sqrt(pow((xx - xc), 2) + pow((yy - yc), 2));

  return (r - rayon(angle(xx, yy), alpha));
}

double Grid::coef(double x, double y)
{
  return 1;
}

double Grid::courbure(double thetha, VectorXd alpha)
{

  double r, r_p, r_s, x_n, p, x_d, rr, rr_p;
  r = rayon(thetha, alpha);
  r_p = rprime(thetha, alpha);
  r_s = rdoubleprime(thetha, alpha);
  p = 1.5;
  rr = pow(r, 2);
  rr_p = pow(r_p, 2);
  x_n = (2 * rr_p + rr - r * r_s);
  x_d = pow((rr_p + rr), p);

  return (x_n) / (x_d);
}

void Grid::BuildT(Eigen::VectorXd alpha)
{

  _tx2.resize(_N[0], _N[1]);
  _ty2.resize(_N[0], _N[1]);
  _xix2_2.resize(_N[0], _N[1]);
  _xiy2_2.resize(_N[0], _N[1]);
  _nbsigne = 0;

  _tx2.setConstant(0);
  _ty2.setConstant(0);
  _xix2_2.setConstant(0);
  _xiy2_2.setConstant(0);

  for (unsigned short j = 0; j < ((unsigned short)(_N[1])); j++)
  {
    for (unsigned short i = 1; i < ((unsigned short)(_N[0] - 1)); i++)
    {

      if (Levelset(i, j, alpha) * Levelset(i + 1, j, alpha) <= 0.)
      {
        _xix2_2(i, j) = 1;
        _nbsigne = _nbsigne + 1;
        _tx2(i, j) = _nbsigne;
      }
      else
      {
        _xix2_2(i, j) = 0;
      }
    }
  }

  for (unsigned short j = 1; j < ((unsigned short)(_N[1] - 1)); j++)
  {
    for (unsigned short i = 0; i < ((unsigned short)(_N[0])); i++)
    {

      if (Levelset(i, j, alpha) * Levelset(i, j + 1, alpha) <= 0.)
      {
        _xiy2_2(i, j) = 1;
        _nbsigne = _nbsigne + 1;
        _ty2(i, j) = _nbsigne;
      }
      else
      {
        _xiy2_2(i, j) = 0;
      }
    }
  }
}

uint Grid::ind(uint i, uint j)
{

  if (i >= _N[0] or j >= _N[1])
  {
    std::cerr << "Grid::index1D: given i,j " << i << "," << j << " but grid has size "
              << _N[0] << "," << _N[1] << std::endl;
    abort();
  }
  return i * _N[1] + (j + 1) + _N[2] + _nbsigne;
}

double Grid::angle_moin(int i)
{
  double pi = acos(-1.);
  double angle_m;

  angle_m = (i)*2. * pi / _N[2] - pi;

  return angle_m;
}

double Grid::angleplus_m(int i)
{
  double pi = acos(-1.);
  double angle_p_m;

  angle_p_m = (i)*2. * pi / _N[2] - pi + 1;

  return angle_p_m;
}

double Grid::Funczone(double theta, Eigen::VectorXd alpha,
                      Eigen::VectorXd angle_debut, Eigen::VectorXd angle_fin)
{
  double theta0, dtheta;
  double funczone;

  double pi = acos(-1.);
  dtheta = pi / (2 * _N[2]);
  funczone = 0;

  for (unsigned short i = 0; i < ((unsigned short)(_N[2])); i++)
  {
    theta0 = angle_debut(i);
    if (angle_fin(i) < theta0)
    {
      if (((theta0 <= theta) && (theta <= pi)) || ((-pi < theta) && (theta <= angle_fin(i))))
        funczone = i + 1;
    }
    else
    {
      if ((theta0 <= theta) && (theta <= angle_fin(i)))
        funczone = i + 1;
    }
  }
  theta0 = 2. * pi - pi;
  if ((theta0 <= theta) && (theta <= (2 * pi - angle_fin(0))))
    funczone = 1;

  return funczone;
}

void Grid::DLevelset(VectorXd alpha)
{
  _u.resize(_N[0], _N[1]);
  _v.resize(_N[0], _N[1]);

  double _norme(0);

  _u.setConstant(0);
  _v.setConstant(0);

  for (unsigned short i = 1; i < ((unsigned short)(_N[0] - 1)); ++i)
  {
    for (unsigned short j = 1; j < ((unsigned short)(_N[1] - 1)); ++j)
    {

      _u(i, j) = (Levelset(i + 1, j, alpha) - Levelset(i - 1, j, alpha)) / (2. * _h[0]);
      _v(i, j) = (Levelset(i, j + 1, alpha) - Levelset(i, j - 1, alpha)) / (2. * _h[1]);
      _norme = sqrt(pow(_u(i, j), 2.) + pow(_v(i, j), 2.));
      _u(i, j) = _u(i, j) / _norme;
      _v(i, j) = _v(i, j) / _norme;
    }
  }

  for (unsigned short j = 1; j < ((unsigned short)(_N[1] - 1)); ++j)
  {
    _u(0, j) = (Levelset(1, j, alpha) - Levelset(0, j, alpha)) / (_h[0]);
    _v(0, j) = (Levelset(0, j + 1, alpha) - Levelset(0, j - 1, alpha)) / (2. * _h[1]);

    _u(_N[0] - 1, j) = (Levelset(_N[0] - 1, j, alpha) - Levelset(_N[0] - 2., j, alpha)) / (_h[0]);
    _v(_N[0] - 1, j) = (Levelset(_N[0] - 1, j + 1, alpha) - Levelset(_N[0] - 1, j - 1, alpha)) / (2. * _h[1]);
  }

  for (unsigned short i = 1; i < ((unsigned short)(_N[0] - 1)); ++i)
  {
    _u(i, 0) = (Levelset(i + 1, 0, alpha) - Levelset(i - 1, 0, alpha)) / (2. * _h[0]);
    _v(i, 0) = (Levelset(i, 1, alpha) - Levelset(i, 0, alpha)) / (_h[1]);

    _u(i, _N[1] - 1) = (Levelset(i + 1, _N[1] - 1, alpha) - Levelset(i - 1, _N[1] - 1, alpha)) / (2. * _h[0]);
    _v(i, _N[1] - 1) = (Levelset(i, _N[1] - 1, alpha) - Levelset(i, _N[1] - 2, alpha)) / (_h[1]);
  }

  _u(0, 0) = (Levelset(1, 0, alpha) - Levelset(0, 0, alpha)) / (_h[0]);
  _v(0, 0) = (Levelset(0, 1, alpha) - Levelset(0, 0, alpha)) / (_h[1]);

  _u(0, _N[1] - 1) = (Levelset(1, _N[1] - 1, alpha) - Levelset(0, _N[1] - 1, alpha)) / (_h[0]);
  _v(0, _N[1] - 1) = (Levelset(0, _N[1] - 1, alpha) - Levelset(0, _N[1] - 2, alpha)) / (_h[1]);

  _u(_N[0] - 1, 0) = (Levelset(_N[0] - 1, 0, alpha) - Levelset(_N[0] - 2, 0, alpha)) / (_h[0]);
  _v(_N[0] - 1, 0) = (Levelset(_N[0] - 1, 1, alpha) - Levelset(_N[0] - 1, 0, alpha)) / (_h[1]);

  _u(_N[0] - 1, _N[1] - 1) = (Levelset(_N[0] - 1, _N[1] - 1, alpha) - Levelset(_N[0] - 2, _N[1] - 1, alpha)) / (_h[0]);
  _v(_N[0] - 1, _N[1] - 1) = (Levelset(_N[0] - 1, _N[1] - 1, alpha) - Levelset(_N[0] - 1, _N[1] - 2, alpha)) / (_h[1]);
}

double Grid::Delta(int i, int j, Eigen::VectorXd alpha)
{
  double term;
  double epsilo(0.0000001);

  _dxphiplus = (Levelset(i + 1, j, alpha) - Levelset(i, j, alpha)) / _h[0];
  _dxphimoins = (Levelset(i, j, alpha) - Levelset(i - 1, j, alpha)) / _h[0];
  _dxphizero = (Levelset(i + 1, j, alpha) - Levelset(i - 1, j, alpha)) / (2. * _h[0]);

  _dyphiplus = (Levelset(i, j + 1, alpha) - Levelset(i, j, alpha)) / _h[1];
  _dyphimoins = (Levelset(i, j, alpha) - Levelset(i, j - 1, alpha)) / _h[1];
  _dyphizero = (Levelset(i, j + 1, alpha) - Levelset(i, j - 1, alpha)) / (2. * _h[1]);

  _absgrad = sqrt(_dxphizero * _dxphizero + _dyphizero * _dyphizero + epsilo);

  _delta = 0.;
  if (Levelset(i, j, alpha) * Levelset(i + 1, j, alpha) <= 0.)
  {
    term = abs(Levelset(i + 1, j, alpha) * _dxphizero) / (_h[0] * _h[0] * abs(_dxphiplus * _absgrad));
    _delta = _delta + term;
  }

  if (Levelset(i, j, alpha) * Levelset(i - 1, j, alpha) < 0.)
  {
    term = abs(Levelset(i - 1, j, alpha) * _dxphizero) / (_h[0] * _h[0] * abs(_dxphimoins * _absgrad));
    _delta = _delta + term;
  }

  if (Levelset(i, j, alpha) * Levelset(i, j + 1, alpha) <= 0.)
  {
    term = abs(Levelset(i, j + 1, alpha) * _dyphizero) / (_h[1] * _h[1] * abs(_dyphiplus * _absgrad));
    _delta = _delta + term;
  }

  if (Levelset(i, j, alpha) * Levelset(i, j - 1, alpha) < 0.)
  {
    term = abs(Levelset(i, j - 1, alpha) * _dyphizero) / (_h[1] * _h[1] * abs(_dyphimoins * _absgrad));
    _delta = _delta + term;
  }

  return _delta;
}

void Grid::inter(VectorXd alpha)
{
  xinter.resize(_N[0], _N[1]);
  yinter.resize(_N[0], _N[1]);

  xinter.setConstant(0);
  yinter.setConstant(0);

  double xx, yy, mid, fin, deb, theta, r1, s;
  double F, FF;
  double yc(0.), xc(0.);

  s = 0.;
  for (unsigned short i = 0; i < ((unsigned short)(_N[0] - 1)); ++i)
  {
    for (unsigned short j = 0; j < ((unsigned short)(_N[1] - 1)); ++j)
    {

      // calcul pour l'axe x
      if (_xix2_2(i, j) == 1)
      {
        // dichotomie
        // initial guess
        deb = _x[i];
        fin = _x[i + 1];

        for (int n = 0; n < 1000; ++n)
        {
          mid = (deb + fin) / 2.;
          yy = _y[j];
          xx = deb;
          r1 = sqrt(pow((yy - yc), 2) + pow((xx - xc), 2));
          theta = asin((yy - yc) / r1);
          F = sqrt(pow((xx - xc), 2) + pow((yy - yc), 2)) - rayon(angle(xx, yy), alpha);
          xx = mid;
          r1 = sqrt(pow((yy - yc), 2)) + pow((xx - xc), 2);
          theta = asin((yy - yc) / r1);
          FF = sqrt(pow((xx - xc), 2) + pow((yy - yc), 2)) - rayon(angle(xx, yy), alpha);

          if (F * FF < 0.)
          {
            fin = mid;
          }
          else
          {
            deb = mid;
          }
        }
        xinter(i, j) = mid;

        if (abs((xinter(i, j) - _x[i]) / _h[0]) <= 0.00000001)
        {
          xinter(i, j) = (float(i) + 1) * _h[0] - squaresize;
        }
        if (abs((xinter(i, j) - _x[i + 1]) / _h[0]) <= 0.00000001)
        {
          xinter(i, j) = (float(i) + 1) * _h[0] - squaresize;
        }
      }

      if (_xiy2_2(i, j) == 1)
      {
        // dichotomie
        // initial guess
        deb = _y[j];
        fin = _y[j + 1];
        for (int n = 0; n < 1000; ++n)
        {
          mid = (deb + fin) / 2.;
          xx = _x[i];
          if (_y[j] >= yc)
          {
            yinter(i, j) = yc + sqrt(pow(1, 2) - pow((_x[i] - xc), 2));
          }
          else
          {
            yinter(i, j) = yc - sqrt(pow(1, 2) - pow((_x[i] - xc), 2));
          }
          yy = deb;
          r1 = sqrt(pow((yy - yc), 2) + pow((xx - xc), 2));
          theta = asin((yy - yc) / r1);
          F = sqrt(pow((xx - xc), 2) + pow((yy - yc), 2)) - rayon(angle(xx, yy), alpha);

          yy = mid;
          r1 = sqrt(pow((yy - yc), 2) + pow((xx - xc), 2));
          theta = asin((yy - yc) / r1);
          FF = sqrt(pow((yy - yc), 2) + pow((xx - xc), 2)) - rayon(angle(xx, yy), alpha);

          if (F * FF <= 0.)
          {
            fin = mid;
          }
          else
          {
            deb = mid;
          }
        }
        yinter(i, j) = mid;

        if (abs((yinter(i, j) - _y[j]) / _h[1]) <= 0.00000001)
        {
          yinter(i, j) = (float(j) + 1) * _h[1] - squaresize;
        }
        if (abs((yinter(i, j) - _y[j + 1]) / _h[1]) <= 0.00000001)
        {
          yinter(i, j) = (float(j) + 1) * _h[1] - squaresize;
        }
      }
    }
  }
}

void Grid::Zonelectrode(Eigen::VectorXd alpha,
                        Eigen::VectorXd angle_debut, Eigen::VectorXd angle_fin)
{

  zone.resize(_nbsigne + 1);
  double r2, theta;
  int xc(0), yc(0);
  double pi = acos(-1.);

  zone.setConstant(0);

  for (unsigned short i = 0; i < ((unsigned short)(_N[0])); i++)
  {
    for (unsigned short j = 0; j < ((unsigned short)(_N[1])); j++)
    {

      if (_tx2(i, j) != 0)
      {
        r2 = sqrt(pow((xinter(i, j) - xc), 2) + pow((_y[j] - yc), 2));
        theta = acos((xinter(i, j) - xc) / r2);
        if (_y[j] <= 0.)
          theta = 2. * pi - theta;

        zone[_tx2(i, j)] = Funczone(angle(xinter(i, j), _y[j]), alpha, angle_debut, angle_fin);
      }

      if (_ty2(i, j) != 0)
      {
        r2 = sqrt(pow((_x[i] - xc), 2) + pow((yinter(i, j) - yc), 2));
        theta = acos((_x[i] - xc) / r2);
        if (yinter(i, j) <= 0.)
          theta = 2. * pi - theta;

        zone[_ty2(i, j)] = Funczone(angle(_x[i], yinter(i, j)), alpha, angle_debut, angle_fin);
      }
    }
  }
}

long double Grid::source(double x, double y)
{
  return -coef(x, y) * (sin(x * y) * (x * x + y * y));
}

long double Grid::val_du(double x, double y, double fix, double fiy)
{
  return coef(x, y) * (y * fix + x * fiy) * cos(x * y);
}

long double Grid::cond_inter(double x, double y)
{
  return sin(x * y);
}

long double Grid::val_dt(double x, double y, Eigen::VectorXd alpha)
{
  return (cos(x * y) / rho(angle(x, y), alpha)) * (cos(angle(x, y)) * y * rprime(angle(x, y), alpha) +
                                                   sin(angle(x, y)) * x * rprime(angle(x, y), alpha) + pow(x, 2) - pow(y, 2));
}

long double Grid::val_dt_norm(double x, double y, Eigen::VectorXd alpha)
{
  return (cos(x * y) / rho(angle(x, y), alpha)) * (x * y + y * x);
}

double Grid::sign_fonc(double x)
{
  double s;

  if (x >= 0)
  {
    s = 1;
  }
  else
  {
    s = -1;
  }
  return s;
}
/////////
void Grid::coord_int()
{

  int i, j;
  coord_inter.resize(2, _nbsigne);

  for (i = 1; i < _N[0] - 1; i++)
  {
    for (j = 1; j < _N[1] - 1; j++)
    {
      if (_tx2(i, j) != 0)
      {
        coord_inter(0, _tx2(i, j) - 1) = xinter(i, j);
        coord_inter(1, _tx2(i, j) - 1) = _y[j];
      }
      if (_ty2(i, j) != 0)
      {
        coord_inter(0, _ty2(i, j) - 1) = _x[i];
        coord_inter(1, _ty2(i, j) - 1) = yinter(i, j);
      }
    }
  }
}
/////////
void Grid::ordre()
{

  uint i, j, k;
  double tmp, tmpp, tmp_0, tmp_2;
  angles_interface.resize(3, _nbsigne);
  indice_interface.resize(_nbsigne);
  indice_final.resize(_nbsigne);
  ind_0.resize(_nbsigne);
  ind_2.resize(_nbsigne);
  x_0.resize(_nbsigne);
  y_0.resize(_nbsigne);
  x_2.resize(_nbsigne);
  y_2.resize(_nbsigne);

  for (i = 1; i < _N[0] - 1; i++)
  {
    for (j = 1; j < _N[1] - 1; j++)
    {
      if (_tx2(i, j) != 0)
      {
        angles_interface(0, _tx2(i, j) - 1) = angle(xinter(i, j), _y[j]);
        indice_interface(_tx2(i, j) - 1) = _tx2(i, j) - 1;
        angles_interface(1, _tx2(i, j) - 1) = xinter(i, j);
        angles_interface(2, _tx2(i, j) - 1) = _y[j];
      }
      if (_ty2(i, j) != 0)
      {
        angles_interface(0, _ty2(i, j) - 1) = angle(_x[i], yinter(i, j));
        indice_interface(_ty2(i, j) - 1) = _ty2(i, j) - 1;
        angles_interface(1, _ty2(i, j) - 1) = _x[i];
        angles_interface(2, _ty2(i, j) - 1) = yinter(i, j);
      }
    }
  }

  for (j = 0; j < _nbsigne - 1; j++)
  {
    for (i = 0; i < _nbsigne - 1; i++)
    {
      if (angles_interface(0, i) > angles_interface(0, i + 1))
      {

        tmp = angles_interface(0, i);
        angles_interface(0, i) = angles_interface(0, i + 1);
        angles_interface(0, i + 1) = tmp;

        tmpp = indice_interface(i);
        indice_interface(i) = indice_interface(i + 1);
        indice_interface(i + 1) = tmpp;

        tmp_0 = angles_interface(1, i);
        angles_interface(1, i) = angles_interface(1, i + 1);
        angles_interface(1, i + 1) = tmp_0;

        tmp_2 = angles_interface(2, i);
        angles_interface(2, i) = angles_interface(2, i + 1);
        angles_interface(2, i + 1) = tmp_2;
      }
    }
  }

  ind_0[0] = indice_interface(_nbsigne - 1);
  ind_2[0] = indice_interface(1);
  x_0[0] = angles_interface(1, _nbsigne - 1);
  y_0[0] = angles_interface(2, _nbsigne - 1);
  x_2[0] = angles_interface(1, 1);
  y_2[0] = angles_interface(2, 1);

  ind_0[_nbsigne - 1] = indice_interface(_nbsigne - 2);
  ind_2[_nbsigne - 1] = indice_interface(0);
  x_0[_nbsigne - 1] = angles_interface(1, _nbsigne - 2);
  y_0[_nbsigne - 1] = angles_interface(2, _nbsigne - 2);
  x_2[_nbsigne - 1] = angles_interface(1, 0);
  y_2[_nbsigne - 1] = angles_interface(2, 0);

  for (k = 1; k < _nbsigne - 1; k++)
  {
    ind_0[k] = indice_interface(k - 1);
    ind_2[k] = indice_interface(k + 1);
    x_0[k] = angles_interface(1, k - 1);
    y_0[k] = angles_interface(2, k - 1);
    x_2[k] = angles_interface(1, k + 1);
    y_2[k] = angles_interface(2, k + 1);
  }

  for (i = 0; i < _nbsigne; i++)
  {
    indice_final(indice_interface(i)) = i;
  }
}

double Grid::longeur_arc(double theta1, double theta2, VectorXd alpha)
{

  double somme_integrale, dt;
  int n, i;

  n = 1000;
  dt = (theta2 - theta1) / n;
  somme_integrale = 0;
  for (i = 0; i < n; i++)
  {
    somme_integrale = somme_integrale + (dt * rho(theta1 + ((i + 1 / 2) * dt), alpha));
  }

  return somme_integrale;
}

void Grid::discret_scheme(uint i, uint j, double fix, double fiy, double x, double y, MatrixXd &M, MatrixXd &M_ind)
{

  int k, l, s_n, s_n1, s_nn, s_n2, IND_1, IND_2;
  double a, aa;
  int som_c1, som_c2, s_indice;
  int carre_1[4], carre_2[4], ind_carre1[4], ind_carre11[4], ind_carre2[4], ind_carre22[4], Tab_ind1[4], Tab_ind2[4];
  double deter1[4], deter2[4];
  MatrixXd norm, det;
  VectorXd a_n;
  int ind_K, ind_I, ind_K2;

  M.resize(6, 2);
  M_ind.resize(6, 2);
  norm.resize(6, 2);
  a_n.resize(6);

  for (k = 0; k < 6; k++)
  {
    norm(k, 0) = 0;
    norm(k, 1) = 0;
    a_n(k) = 0;
  }

  for (k = 0; k < 6; k++)
  {
    norm(k, 0) = (x - M(k, 0)) / sqrt(pow((x - M(k, 0)), 2) + pow((y - M(k, 1)), 2));
    norm(k, 1) = (y - M(k, 1)) / sqrt(pow((x - M(k, 0)), 2) + pow((y - M(k, 1)), 2));
  }
  // produit scalaire avec les composante normales
  for (k = 0; k < 6; k++)
  {
    a_n(k) = (norm(k, 0) * fix + norm(k, 1) * fiy);
  }
  // initialise les sommes
  som_c1 = 0;
  som_c2 = 0;
  for (l = 0; l < 4; l++)
  {
    carre_1[l] = 0;
    carre_2[l] = 0;
  }
  // remplir les deux carré
  for (l = 0; l < 4; l++)
  {
    if (a_n(l) > 0)
    {
      carre_1[l] = 1;
    }
    else
    {
      carre_1[l] = 0;
    }
    som_c1 = som_c1 + carre_1[l];
  }
  // remplir le deuxième carré
  for (l = 2; l < 6; l++)
  {
    if (a_n(l) > 0)
    {
      carre_2[l - 2] = 1;
    }
    else
    {
      carre_2[l - 2] = 0;
    }
    som_c2 = som_c2 + carre_2[l - 2];
  }

  for (l = 0; l < 4; l++)
  {
    ind_carre1[l] = 0;
  }

  if (som_c1 > som_c2)
  {

    s_indice = som_c1;
    s_n = 0;
    for (l = 0; l < 4; l++)
    {
      if (carre_1[l] != 0)
      {
        ind_carre1[s_n] = l;
        s_n = s_n + 1;
      }
    }

    // 2ieme etape du travaille dans le cas ou carre1 gagne la somme
    // deter1.resize(s_n);
    if (s_indice != 2)
    {
      // je choisit xk,yk
      for (k = 0; k < 4; k++)
      {
        deter1[k] = 0;
      }
      a = 0;
      for (k = 0; k < s_n; k++)
      {
        deter1[k] = norm(ind_carre1[k], 0) * fiy - norm(ind_carre1[k], 1) * fix;
        a = a + sign_fonc(deter1[k]);
      }

      if (a == 1 || a == -1.)
      {
        for (k = 0; k < s_n; k++)
        {
          if (a == 1.)
          {
            if (sign_fonc(deter1[k]) == -1.)
            {
              x_k = M(ind_carre1[k], 0);
              y_k = M(ind_carre1[k], 1);
              ind_K = ind_carre1[k];
              indice_k = ind(M_ind(ind_K, 0), M_ind(ind_K, 1));
            }
          }
          else if (a == -1)
          {
            if (sign_fonc(deter1[k]) == 1.)
            {
              x_k = M(ind_carre1[k], 0);
              y_k = M(ind_carre1[k], 1);
              ind_K = ind_carre1[k];
              indice_k = ind(M_ind(ind_K, 0), M_ind(ind_K, 1));
            }
          }
        }
      }
      else
      {
        cout << "Im not theoretically well: case a=" << a << endl;
        x_k = M(ind_carre1[0], 0);
        y_k = M(ind_carre1[0], 1);
        x_i = M(ind_carre1[1], 0);
        y_i = M(ind_carre1[1], 1);
        indice_k = ind(M_ind(ind_carre1[0], 0), M_ind(ind_carre1[0], 1));
        indice_k = ind(M_ind(ind_carre1[1], 0), M_ind(ind_carre1[1], 1));
      }

      // je choisit xi,yi
      for (k = 0; k < 4; k++)
      {
        ind_carre11[k] = 0;
      }
      s_n1 = 0;
      for (k = 0; k < s_n; k++)
      {
        if (ind_carre1[k] != ind_K)
        {
          ind_carre11[s_n1] = ind_carre1[k];
          s_n1 = s_n1 + 1;
        }
      }

      for (k = 0; k < s_n1 - 1; k++)
      {
        if (a_n(ind_carre11[k]) > a_n(ind_carre11[k + 1]))
        {
          x_i = M(ind_carre11[k], 0);
          y_i = M(ind_carre11[k], 1);
          indice_i = ind(M_ind(ind_carre11[k], 0), M_ind(ind_carre11[k], 1));
        }
        else
        {
          x_i = M(ind_carre11[k + 1], 0);
          y_i = M(ind_carre11[k + 1], 1);
          indice_i = ind(M_ind(ind_carre11[k + 1], 0), M_ind(ind_carre11[k + 1], 1));
        }
      }
    }

    if (s_indice == 2)
    { // cas ou som_c1=2 et somc_2=1
      for (k = 0; k < 1; k++)
      {
        x_k = M(ind_carre1[k], 0);
        y_k = M(ind_carre1[k], 1);
        x_i = M(ind_carre1[k + 1], 0);
        y_i = M(ind_carre1[k + 1], 1);
        indice_k = ind(M_ind(ind_carre1[k], 0), M_ind(ind_carre1[k], 1));
        ;
        indice_i = ind(M_ind(ind_carre1[k + 1], 0), M_ind(ind_carre1[k + 1], 1));
      }
    }
  }

  for (l = 0; l < 4; l++)
  {
    ind_carre2[l] = 0;
  }
  if (som_c2 > som_c1)
  {
    s_indice = som_c2;
    s_nn = 0;
    for (l = 0; l < 4; l++)
    {
      if (carre_2[l] != 0)
      {
        ind_carre2[s_nn] = l + 2;
        s_nn = s_nn + 1;
      }
    }

    // 2ieme etape du travaille dans le cas ou carre1 gagne la somme

    if (s_indice != 2)
    {
      for (k = 0; k < s_nn; k++)
      {
        deter2[k] = 0;
      }

      // je choisit xk,yk
      aa = 0;
      for (k = 0; k < s_nn; k++)
      {
        deter2[k] = norm(ind_carre2[k], 0) * fiy - norm(ind_carre2[k], 1) * fix;
        aa = aa + sign_fonc(deter2[k]);
      }

      if (aa == 1 || aa == -1)
      {
        for (k = 0; k < s_nn; k++)
        {
          if (aa == 1)
          {
            if (sign_fonc(deter2[k]) == -1)
            {
              x_k = M(ind_carre2[k], 0);
              y_k = M(ind_carre2[k], 1);
              ind_K2 = ind_carre2[k];
              indice_k = ind(M_ind(ind_K2, 0), M_ind(ind_K2, 1));
            }
          }
          else if (aa == -1)
          {
            if (sign_fonc(deter2[k]) == 1)
            {
              x_k = M(ind_carre2[k], 0);
              y_k = M(ind_carre2[k], 1);
              ind_K2 = ind_carre2[k];
              indice_k = ind(M_ind(ind_K2, 0), M_ind(ind_K2, 1));
            }
          }
        }
      }
      else
      {
        cout << "Im not theoretically well: case aa=" << aa << endl;
        x_k = M(ind_carre2[0], 0);
        y_k = M(ind_carre2[0], 1);
        x_i = M(ind_carre2[1], 0);
        y_i = M(ind_carre2[1], 1);
        indice_k = ind(M_ind(ind_carre2[0], 0), M_ind(ind_carre2[0], 1));
        indice_k = ind(M_ind(ind_carre2[1], 0), M_ind(ind_carre2[1], 1));
      }

      // je choisit xi,yi
      for (k = 0; k < 4; k++)
      {
        ind_carre22[k] = 0;
      }
      s_n2 = 0;
      for (k = 0; k < s_nn; k++)
      {
        if (ind_carre2[k] != ind_K2)
        {
          ind_carre22[s_n2] = ind_carre2[k];
          s_n2 = s_n2 + 1;
        }
      }

      for (k = 0; k < s_n2 - 1; k++)
      {
        if (a_n(ind_carre22[k]) > a_n(ind_carre22[k + 1]))
        {
          x_i = M(ind_carre22[k], 0);
          y_i = M(ind_carre22[k], 1);
          indice_i = ind(M_ind(ind_carre22[k], 0), M_ind(ind_carre22[k], 1));
        }
        else
        {
          x_i = M(ind_carre22[k + 1], 0);
          y_i = M(ind_carre22[k + 1], 1);
          indice_i = ind(M_ind(ind_carre22[k + 1], 0), M_ind(ind_carre22[k + 1], 1));
        }
      }
    }
    if (s_indice == 2)
    { // cas ou som_c1=2 et somc_2=1
      for (k = 0; k < 1; k++)
      {
        x_k = M(ind_carre2[k], 0);
        y_k = M(ind_carre2[k], 1);
        x_i = M(ind_carre2[k + 1], 0);
        y_i = M(ind_carre2[k + 1], 1);
        indice_k = ind(M_ind(ind_carre2[k], 0), M_ind(ind_carre2[k], 1));
        indice_i = ind(M_ind(ind_carre2[k + 1], 0), M_ind(ind_carre2[k + 1], 1));
      }
    }
  }

  for (l = 0; l < 4; l++)
  {
    Tab_ind1[l] = 0;
    Tab_ind2[l] = 0;
  }
  if (som_c1 == som_c2)
  {

    IND_1 = 0;
    for (l = 0; l < 4; l++)
    {
      if (carre_1[l] != 0)
      {
        Tab_ind1[IND_1] = l;
        IND_1 = IND_1 + 1;
      }
    }

    IND_2 = 0;
    for (l = 0; l < 4; l++)
    {
      if (carre_1[l] != 0)
      {
        Tab_ind2[IND_2] = l + 2;
        IND_2 = IND_2 + 1;
      }
    }

    det.resize(IND_1, 2);
    for (k = 0; k < IND_1; k++)
    {
      det(k, 0) = norm(Tab_ind1[k], 0) * fiy - norm(Tab_ind1[k], 1) * fix;
      det(k, 1) = norm(Tab_ind2[k], 0) * fiy - norm(Tab_ind2[k], 1) * fix;
    }

    if (det(0, 0) * det(1, 0) < 0)
    {
      for (k = 0; k < IND_1 - 1; k++)
      {
        x_k = M(Tab_ind1[k], 0);
        y_k = M(Tab_ind1[k], 1);
        x_i = M(Tab_ind1[k + 1], 0);
        y_i = M(Tab_ind1[k + 1], 1);
        indice_k = ind(M_ind(Tab_ind1[k], 0), M_ind(Tab_ind1[k], 1));
        indice_i = ind(M_ind(Tab_ind1[k + 1], 0), M_ind(Tab_ind1[k + 1], 1));
      }
    }

    if (det(0, 1) * det(1, 1) < 0)
    {
      for (k = 0; k < IND_2 - 1; k++)
      {
        x_k = M(Tab_ind2[k], 0);
        y_k = M(Tab_ind2[k], 1);
        x_i = M(Tab_ind2[k + 1], 0);
        y_i = M(Tab_ind2[k + 1], 1);
        indice_k = ind(M_ind(Tab_ind2[k], 0), M_ind(Tab_ind2[k], 1));
        indice_i = ind(M_ind(Tab_ind2[k + 1], 0), M_ind(Tab_ind2[k + 1], 1));
      }
    }
  }
}
