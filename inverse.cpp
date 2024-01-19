#include "inverse.hpp"
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include "Eigen/Eigen/Eigen"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Core"
#include "Eigen/Eigen/Sparse"
#include "Eigen/Eigen/IterativeLinearSolvers"
#include "Eigen/unsupported/Eigen/src/IterativeSolvers/GMRES.h"

using namespace std;

Inverser::Inverser(Parameters &p, HeatProblem &pb)
{
  _params = &p;

  if (_params->hasKey("Nx_coarse"))
    _params->set("Nx", _params->getString("Nx_coarse"));
  if (_params->hasKey("Ny_coarse"))
    _params->set("Ny", _params->getString("Ny_coarse"));

  _pb = new HeatProblem(p);
  _grid = Grid(p);
  _assbl = EITAssembler(&_grid, _pb); // not sure about _pb
  _writer = Writer(p, _grid);

  squaresize = 2.;

  _x.resize(_grid.getNbPoints(0)); // (xmin+h, ..., xmax-h)
  for (int i = 0; i < _grid.getNbPoints(0); ++i)
    _x[i] = ((i + 1) - 0.5) * _grid.getSpacing(0) - squaresize;

  _y.resize(_grid.getNbPoints(1));
  for (int i = 0; i < _grid.getNbPoints(1); ++i)
    _y[i] = ((i + 1) - 0.5) * _grid.getSpacing(1) - squaresize;
}
Inverser::~Inverser()
{
  delete _pb;
}
////////////////////////////////////////////////////////////////////////////////

double Inverser::F_TV(Eigen::VectorXd alpha,
                      MatrixXd Immm,
                      MatrixXd Ummm, double t,
                      Eigen::VectorXd sig,
                      Eigen::VectorXd sig0,
                      Eigen::VectorXd direction,
                      Eigen::VectorXd &angle_debut,
                      Eigen::VectorXd &angle_fin,
                      VectorXd &Z)
{

  double func;
  double uum, uup, xm, xp, ym, yp;
  _grid.BuildT(alpha);
  _grid.inter(alpha);
  uint n = _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2);
  int ne = _grid.getNbPoints(2);

  VectorXd I_mmm(ne);
  double epsilon = 0.001;
  Immm.resize(ne, ne);
  Ummm.resize(ne, ne);
  MatrixXf uuu(ne, n);
  VectorXd gradsigma(n);

  sig.resize(n);
  sig0.resize(n);
  direction.resize(n);

  for (int i = 0; i < _grid.getNbPoints(0); i++)
  {
    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      sig(_grid.ind(i, j) - 1) += t * direction(_grid.ind(i, j) - 1);
    }
  }

  for (int l = 0; l < ne; l++)
  {

    for (int i = 0; i < ne; i++)
    {
      I_mmm(i) = Immm(l, i);
    }
    _pb->run(sig, alpha, I_mmm, angle_debut, angle_fin, Z);
    for (int i = 0; i < n; i++)
    {
      uuu(l, i) = _pb->GetSolution()[i];
    }
  }

  gradsigma.setConstant(0);

  for (int i = 0; i < _grid.getNbPoints(0); i++)
  {
    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      if (_grid.Levelset(i, j, alpha) <= 0)
      {

        uup = (sig(_grid.ind(i + 1, j) - 1) - sig0(_grid.ind(i + 1, j) - 1));
        xp = _x[i + 1];
        if (_grid.get_tx2(i, j) != 0)
        {
          uup = (sig(_grid.get_tx2(i, j) - 1) - sig0(_grid.get_tx2(i, j) - 1));
          xp = _grid.get_xinter(i, j);
        }

        uum = (sig(_grid.ind(i - 1, j) - 1) - sig0(_grid.ind(i - 1, j) - 1));
        xm = _x[i - 1];
        if (_grid.get_tx2(i - 1, j) != 0)
        {
          uum = (sig(_grid.get_tx2(i - 1, j) - 1) - sig0(_grid.get_tx2(i - 1, j) - 1));
          xm = _grid.get_xinter(i - 1, j);
        }
        gradsigma(_grid.ind(i, j) - 1) += ((uup - uum) / (xp - xm));

        uup = (sig(_grid.ind(i, j + 1) - 1) - sig0(_grid.ind(i, j + 1) - 1));
        yp = _y[j + 1];
        if (_grid.get_ty2(i, j) != 0)
        {
          uup = (sig(_grid.get_ty2(i, j) - 1) - sig0(_grid.get_ty2(i, j) - 1));
          yp = _grid.get_yinter(i, j);
        }

        uum = (sig(_grid.ind(i, j - 1) - 1) - sig0(_grid.ind(i, j - 1) - 1));
        ym = _y[j - 1];
        if (_grid.get_ty2(i, j - 1) != 0)
        {
          uum = (sig(_grid.get_ty2(i, j - 1) - 1) - sig0(_grid.get_ty2(i, j - 1) - 1));
          ym = _grid.get_yinter(i, j - 1);
        }
        gradsigma(_grid.ind(i, j) - 1) += ((uup - uum) / (yp - ym));
      }
    }
  }

  double nom = gradsigma.cwiseAbs().dot(gradsigma.cwiseAbs());
  double den = sqrt(gradsigma.cwiseAbs().dot(gradsigma.cwiseAbs()) + 0.001);

  func = 0;

  for (int l = 0; l < ne; l++)
  {
    for (int m = 0; m < ne; m++)
    {
      func = func + 0.5 * pow((Ummm(l, m) - uuu(l, _grid.Get_nbsigne() + m)), 2);
    }
  }

  for (int i = 0; i < _grid.getNbPoints(0); i++)
  {
    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      if (_grid.Levelset(i, j, alpha) <= 0)
      {
        // add the integral
        func += epsilon / 2 * (nom / den) * _grid.getSpacing(0) * _grid.getSpacing(1);
      }
    }
  }

  return func;
}

////////////////////////////////////////////////////////////////////////////////

double Inverser::F(Eigen::VectorXd alpha,
                   MatrixXd Immm,
                   MatrixXd Ummm, double t,
                   Eigen::VectorXd sig,
                   Eigen::VectorXd sig0,
                   Eigen::VectorXd direction,
                   Eigen::VectorXd &angle_debut,
                   Eigen::VectorXd &angle_fin,
                   VectorXd &Z)
{

  double func;
  double uum, uup, xm, xp, ym, yp;
  _grid.BuildT(alpha);
  _grid.inter(alpha);
  uint n = _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2);
  int ne = _grid.getNbPoints(2);

  VectorXd I_mmm(ne);
  double epsilon = 0.001;
  Immm.resize(ne, ne);
  Ummm.resize(ne, ne);
  MatrixXf uuu(ne, n);
  MatrixXf grads(_grid.getNbPoints(0), _grid.getNbPoints(1));

  sig.resize(n);
  sig0.resize(n);
  direction.resize(n);

  for (int i = 0; i < _grid.getNbPoints(0); i++)
  {
    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      sig(_grid.ind(i, j) - 1) += t * direction(_grid.ind(i, j) - 1);
    }
  }

  for (int l = 0; l < ne; l++)
  {

    for (int i = 0; i < ne; i++)
    {
      I_mmm(i) = Immm(l, i);
    }
    _pb->run(sig, alpha, I_mmm, angle_debut, angle_fin, Z);
    for (int i = 0; i < n; i++)
    {
      uuu(l, i) = _pb->GetSolution()[i];
    }
  }

  for (int i = 0; i < _grid.getNbPoints(0); i++)
  {
    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      grads(i, j) = 0;
    }
  }

  for (int i = 0; i < _grid.getNbPoints(0); i++)
  {
    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      if (_grid.Levelset(i, j, alpha) <= 0)
      {

        uup = (sig(_grid.ind(i + 1, j) - 1) - sig0(_grid.ind(i + 1, j) - 1));
        xp = _x[i + 1];
        if (_grid.get_tx2(i, j) != 0)
        {
          uup = (sig(_grid.get_tx2(i, j) - 1) - sig0(_grid.get_tx2(i, j) - 1));
          xp = _grid.get_xinter(i, j);
        }

        uum = (sig(_grid.ind(i - 1, j) - 1) - sig0(_grid.ind(i - 1, j) - 1));
        xm = _x[i - 1];
        if (_grid.get_tx2(i - 1, j) != 0)
        {
          uum = (sig(_grid.get_tx2(i - 1, j) - 1) - sig0(_grid.get_tx2(i - 1, j) - 1));
          xm = _grid.get_xinter(i - 1, j);
        }
        grads(i, j) = grads(i, j) + ((uup - uum) / (xp - xm));

        uup = (sig(_grid.ind(i, j + 1) - 1) - sig0(_grid.ind(i, j + 1) - 1));
        yp = _y[j + 1];
        if (_grid.get_ty2(i, j) != 0)
        {
          uup = (sig(_grid.get_ty2(i, j) - 1) - sig0(_grid.get_ty2(i, j) - 1));
          yp = _grid.get_yinter(i, j);
        }

        uum = (sig(_grid.ind(i, j - 1) - 1) - sig0(_grid.ind(i, j - 1) - 1));
        ym = _y[j - 1];
        if (_grid.get_ty2(i, j - 1) != 0)
        {
          uum = (sig(_grid.get_ty2(i, j - 1) - 1) - sig0(_grid.get_ty2(i, j - 1) - 1));
          ym = _grid.get_yinter(i, j - 1);
        }
        grads(i, j) = grads(i, j) + ((uup - uum) / (yp - ym));
      }
    }
  }

  func = 0;

  for (int l = 0; l < ne; l++)
  {
    for (int m = 0; m < ne; m++)
    {
      func = func + 0.5 * pow((Ummm(l, m) - uuu(l, _grid.Get_nbsigne() + m)), 2);
    }
  }

  for (int i = 0; i < _grid.getNbPoints(0); i++)
  {
    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      if (_grid.Levelset(i, j, alpha) <= 0)
      {
        func += _grid.getSpacing(0) * _grid.getSpacing(1) * 0.5 *
                epsilon * (pow(sig(_grid.ind(i, j) - 1) - sig0(_grid.ind(i, j) - 1), 2) + pow(grads(i, j), 2));
      }
    }
  }

  return func;
}

double Inverser::cost_function(VectorXd &alpha,
                               VectorXd &angle_debut,
                               VectorXd &angle_fin,
                               VectorXd &sigma,
                               VectorXd &contact,
                               VectorXd &contact_terre,
                               MatrixXd &current_patern,
                               MatrixXd &measured_potentials,
                               double t,
                               VectorXd &dir_contact)
{
  double cost = 0;

  uint n = _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2);
  int ne = _grid.getNbPoints(2);

  VectorXd Im(ne);
  double epsilon = 0.001;
  MatrixXf eit_sol(ne, n);

  // sigma.resize(n);
  // sigma.setConstant(1);

  contact -= t * dir_contact;

  for (int l = 0; l < ne; l++)
  {

    for (int i = 0; i < ne; i++)
    {
      Im(i) = current_patern(l, i);
    }
    _pb->run(sigma, alpha, Im, angle_debut, angle_fin, contact);
    for (int i = 0; i < n; i++)
    {
      eit_sol(l, i) = _pb->GetSolution()[i];
    }
  }

  cost = 0;

  for (int l = 0; l < ne - 1; l++)
  {
    for (int m = 0; m < ne; m++)
    {

      double U = (measured_potentials(l, m) - eit_sol(l, _grid.Get_nbsigne() + m));
      double Z = contact(m) - contact_terre(m);
      cost += 0.5 * U * U + (epsilon / 2) * Z * Z;
    }
  }

  contact += t * dir_contact;
  return cost;
}

double Inverser::cost_function_both(VectorXd &alpha,
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
                                    VectorXd &dir_conductivity)

{
  double cost = 0;

  uint n = _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2);
  int ne = _grid.getNbPoints(2);

  VectorXd Im(ne);
  double epsilon = 0.001;
  MatrixXf eit_sol(ne, n);
  double uum, uup, xm, xp, ym, yp;
  MatrixXf grads(_grid.getNbPoints(0), _grid.getNbPoints(1));

  sigma.resize(n);
  sigma.setConstant(1);

  contact -= t * dir_contact;
  sigma -= t * dir_conductivity;

  for (int l = 0; l < ne; l++)
  {

    for (int i = 0; i < ne; i++)
    {
      Im(i) = current_patern(l, i);
    }
    _pb->run(sigma, alpha, Im, angle_debut, angle_fin, contact);
    for (int i = 0; i < n; i++)
    {
      eit_sol(l, i) = _pb->GetSolution()[i];
    }
  }

  grads.setConstant(0);
  for (int i = 0; i < _grid.getNbPoints(0); i++)
  {
    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      if (_grid.Levelset(i, j, alpha) <= 0)
      {

        uup = (sigma(_grid.ind(i + 1, j) - 1) - sigma0(_grid.ind(i + 1, j) - 1));
        xp = _x[i + 1];
        if (_grid.get_tx2(i, j) != 0)
        {
          uup = (sigma(_grid.get_tx2(i, j) - 1) - sigma0(_grid.get_tx2(i, j) - 1));
          xp = _grid.get_xinter(i, j);
        }

        uum = (sigma(_grid.ind(i - 1, j) - 1) - sigma0(_grid.ind(i - 1, j) - 1));
        xm = _x[i - 1];
        if (_grid.get_tx2(i - 1, j) != 0)
        {
          uum = (sigma(_grid.get_tx2(i - 1, j) - 1) - sigma0(_grid.get_tx2(i - 1, j) - 1));
          xm = _grid.get_xinter(i - 1, j);
        }
        grads(i, j) = grads(i, j) + ((uup - uum) / (xp - xm));

        uup = (sigma(_grid.ind(i, j + 1) - 1) - sigma0(_grid.ind(i, j + 1) - 1));
        yp = _y[j + 1];
        if (_grid.get_ty2(i, j) != 0)
        {
          uup = (sigma(_grid.get_ty2(i, j) - 1) - sigma0(_grid.get_ty2(i, j) - 1));
          yp = _grid.get_yinter(i, j);
        }

        uum = (sigma(_grid.ind(i, j - 1) - 1) - sigma0(_grid.ind(i, j - 1) - 1));
        ym = _y[j - 1];
        if (_grid.get_ty2(i, j - 1) != 0)
        {
          uum = (sigma(_grid.get_ty2(i, j - 1) - 1) - sigma0(_grid.get_ty2(i, j - 1) - 1));
          ym = _grid.get_yinter(i, j - 1);
        }
        grads(i, j) = grads(i, j) + ((uup - uum) / (yp - ym));
      }
    }
  }

  // func = 0;

  // for (int l = 0; l < ne; l++)
  // {
  //   for (int m = 0; m < ne; m++)
  //   {
  //     func = func + 0.5 * pow((Ummm(l, m) - uuu(l, _grid.Get_nbsigne() + m)), 2);
  //   }
  // }

  // for (int i = 0; i < _grid.getNbPoints(0); i++)
  // {
  //   for (int j = 0; j < _grid.getNbPoints(1); j++)
  //   {
  //     if (_grid.Levelset(i, j, alpha) <= 0)
  //     {
  //       func += _grid.getSpacing(0) * _grid.getSpacing(1) * 0.5 *
  //               epsilon * (pow(sig(_grid.ind(i, j) - 1) - sig0(_grid.ind(i, j) - 1), 2) + pow(grads(i, j), 2));
  //     }
  //   }
  // }

  cost = 0;

  for (int l = 0; l < ne - 1; l++)
  {
    for (int m = 0; m < ne; m++)
    {

      double U = (measured_potentials(l, m) - eit_sol(l, _grid.Get_nbsigne() + m));
      double Z = contact(m) - contact_terre(m);

      cost += 0.5 * U * U + (epsilon / 2) * Z * Z;
      ;
    }
  }

  for (int i = 0; i < _grid.getNbPoints(0); i++)
    for (int j = 0; j < _grid.getNbPoints(1); j++)
      if (_grid.Levelset(i, j, alpha) <= 0)
        cost += _grid.getSpacing(0) * _grid.getSpacing(1) * 0.5 * epsilon * (pow(sigma(_grid.ind(i, j) - 1) - sigma0(_grid.ind(i, j) - 1), 2) + pow(grads(i, j), 2));

  return cost;
}

void Inverser::runinvcontact(VectorXd &alpha,
                             VectorXd &angle_debut,
                             VectorXd &angle_fin,
                             VectorXd &sigma)
{
  _grid.BuildT(alpha);
  _grid.inter(alpha);
  _grid.Zonelectrode(alpha, angle_debut, angle_fin);

  uint n = _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2);
  int ne = _grid.getNbPoints(2);

  std::cout << "Nx = " << _grid.getNbPoints(0)
            << ", Ny = " << _grid.getNbPoints(1)
            << std::endl;

  _sigma.resize(n);
  double epsilon = 0.001;

  VectorXd Im(ne);
  VectorXd descent_contact(ne), contact(ne), contact_terre(ne);
  MatrixXd descent(ne - 1, ne);

  descent.setConstant(0);
  descent_contact.setConstant(0);
  contact_terre.setConstant(2);
  contact = contact_terre;

  MatrixXd currents(ne, ne);
  MatrixXd measurments(ne, ne);
  MatrixXd U(ne, ne);
  MatrixXd u(ne, n);
  MatrixXd w(ne, n);

  _sigma.setConstant(1);

  // Ouvrir les fichier et faire la lecture
  std::string folder = _params->getString("output dir");
  for (int l = 0; l < ne - 1; l++)
  {
    std::string nomfichier = folder + "/measure-elec-" + std::to_string(l) + ".txt";
    std::ifstream fichier(nomfichier.c_str()); // le constructeur de std::ifstream n'accepte pas de string
    while (!fichier.eof())
    {
      for (int i = 0; i < ne; i++)
      {
        fichier >> setprecision(17) >> currents(l, i) >> measurments(l, i);
      }
    }

    for (int i = 0; i < ne; i++)
    {
      measurments(l, i) = generator(measurments(l, i));
    }
  }

  for (int t = 0; t < 300; t++)
  {

    descent_contact.setConstant(0);

    for (int i = 0; i < ne; i++)
    {
      descent_contact(i) = epsilon * (contact(i) - contact_terre(i));
    }

    for (int l = 0; l < ne - 1; l++)
    {

      for (int i = 0; i < ne; i++)
      {
        Im(i) = currents(l, i);
      }
      _pb->run(_sigma, alpha, Im, angle_debut, angle_fin, contact);
      for (int i = 0; i < n; i++)
      {
        u(l, i) = _pb->GetSolution()[i];
      }

      for (int i = 0; i < ne; i++)
      {
        U(l, i) = u(l, _grid.Get_nbsigne() + i);
        Im(i) = U(l, i) - measurments(l, i);
      }

      _pb->run(_sigma, alpha, Im, angle_debut, angle_fin, contact);
      for (int i = 0; i < n; i++)
      {
        w(l, i) = _pb->GetSolution()[i];
      }

      for (int i = 0; i < ne; i++)
      {
        descent(l, i) = _assbl.Compute_integral(u, w, alpha, angle_debut, angle_fin, i, l);
      }
    }

    for (int i = 0; i < ne; i++)
    {
      for (int l = 0; l < ne - 1; l++)
      {
        descent_contact(i) -= descent(l, i);
      }
    }
    // normalisation
    double norm = 0;
    for (int i = 0; i < ne; i++)
      norm += descent_contact(i) * descent_contact(i);

    norm = sqrt(norm);

    descent_contact = descent_contact / norm;

    descent_contact = 0.05 * descent_contact;

    for (int i = 0; i < ne; i++)
    {
      contact(i) = contact(i) - descent_contact(i);
      std::cout << contact(i) << "  " << descent_contact(i) << "  " << i << std::endl;
    }

    double F = cost_function(alpha,
                             angle_debut,
                             angle_fin,
                             sigma,
                             contact,
                             contact_terre,
                             currents,
                             measurments,
                             0.05,
                             descent_contact);

    std::cout << "Fcontact::"
              << "  " << F << " " << t << std::endl;
  }
}

double Inverser::generator(float vect_meas)
{ /* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */
  double temp1;
  double temp2;
  double result;
  int p;
  double pi = 3.1415926536;

  p = 1;

  while (p > 0)
  {
    temp2 = (rand() / ((double)RAND_MAX));
    /*  rand() function generates an
        integer between 0 and  RAND_MAX,
         which is defined in stdlib.h.
    */

    if (temp2 == 0)
    { // temp2 is >= (RAND_MAX / 2)
      p = 1;
    } // end if
    else
    { // temp2 is < (RAND_MAX / 2)
      p = -1;
    } // end else

  } // end while()

  temp1 = cos((2.0 * pi) * rand() / ((double)RAND_MAX));
  result = sqrt(-2.0 * log(temp2)) * temp1;

  double noise = _params->getReal("noise");
  return (vect_meas + (noise * result)); // return the generated random sample to the caller

} // end AWGN_generator()

////////////////////////////////////////////////////////////////////////////////////

void Inverser::runinv(VectorXd &alpha,
                      VectorXd &angle_debut,
                      VectorXd &angle_fin,
                      VectorXd &Z)
{
  _grid.BuildT(alpha);
  _grid.inter(alpha);
  uint n = _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2);
  int ne = _grid.getNbPoints(2);
  std::cout << "Nx = " << _grid.getNbPoints(0) << ", Ny = " << _grid.getNbPoints(1) << std::endl;
  Eigen::VectorXd Im, qqq;
  MatrixXd aa(ne, ne);
  MatrixXd bb(ne, ne);
  MatrixXd Imm(ne, ne);
  MatrixXd Umm(ne, ne);
  MatrixXd Umm2(ne, ne);
  MatrixXd u(ne, n);
  MatrixXd w(ne, n);
  MatrixXf dsig(ne, n);
  MatrixXf grads(_grid.getNbPoints(0), _grid.getNbPoints(1));
  MatrixXf ppp(_grid.getNbPoints(0), _grid.getNbPoints(1));
  double sumu, uum, uup, xm, xp, ym, yp;
  double t_moin, t_plus, t_l, t_r, tt, a, fr, fl;

  // std::cout<<std::setprecision(17)<< _grid.Get_nbsigne()<<std::endl;

  _sigma.resize(n);
  _sigma0.resize(n);
  sources.resize(n);
  _dir.resize(n);
  qqq.resize(n);

  Im.resize(ne);
  double pi = acos(-1.), Fsigma;
  string ligne, maligne;
  double epsilon = 0.001;

  for (int i = 0; i < n; i++)
  {
    _sigma0[i] = 1;
    _sigma[i] = _sigma0[i];
  }

  // Ouvrir les fichier et faire la lecture
  std::string folder = _params->getString("output dir");
  for (int l = 0; l < ne - 1; l++)
  {
    std::string nomfichier = folder + "/measure-elec-" + std::to_string(l) + ".txt";
    std::ifstream fichier(nomfichier.c_str()); // le constructeur de std::ifstream n'accepte pas de string
    while (!fichier.eof())
    {
      for (int i = 0; i < ne; i++)
      {
        fichier >> setprecision(17) >> aa(l, i) >> bb(l, i);
      }
    }

    for (int i = 0; i < ne; i++)
    {
      Imm(l, i) = aa(l, i);
      Umm(l, i) = generator(bb(l, i));
    }
  }

  for (int t = 0; t < 20; t++)
  {
    //////////////////////////////////////////////////////////////////////////////
    for (int l = 0; l < ne; l++)
    {

      for (int i = 0; i < ne; i++)
      {
        Im(i) = Imm(l, i);
      }
      _pb->run(_sigma, alpha, Im, angle_debut, angle_fin, Z);
      for (int i = 0; i < n; i++)
      {
        u(l, i) = _pb->GetSolution()[i];
      }

      for (int i = 0; i < ne; i++)
      {
        Umm2(l, i) = u(l, _grid.Get_nbsigne() + i);
        Im(i) = Umm2(l, i) - Umm(l, i);
      }

      _pb->run(_sigma, alpha, Im, angle_debut, angle_fin, Z);
      for (int i = 0; i < n; i++)
      {
        w(l, i) = _pb->GetSolution()[i];
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Builduing the matrix to get the directory
      _assbl.BuildDirMatrix(alpha, angle_debut, angle_fin, _sigma);
      // Builduing the source term
      _assbl.Buildsourcedir1(l, u, w, alpha);
      // _assbl.Buildsourcedir2(_sigma, _sigma0, alpha);

      for (int i = 0; i < _grid.getNbPoints(0); i++)
      {
        for (int j = 0; j < _grid.getNbPoints(1); j++)
        {
          dsig(l, _grid.ind(i, j) - 1) = _assbl.Getdsigma(i, j);
        }
      }
    }

    VectorXd G = _assbl.Build_dirVT_G(alpha, _sigma);
    VectorXd gradG = _assbl.grad_of_vec(alpha, G);
    std::cout << G.size() << "  " << sources.size() << "  " << _sigma.size() << "  " << gradG.size() << std::endl;

    sources.setConstant(0);

    // for (int j = 0; j < _grid.getNbPoints(1); j++)
    // {
    //   for (int i = 0; i < _grid.getNbPoints(0); i++)
    //   {

    //     for (int l = 0; l < ne; l++)
    //     {
    //       sources(_grid.ind(i, j) - 1) = sources(_grid.ind(i, j) - 1) - dsig(l, _grid.ind(i, j) - 1);
    //     }
    //     sources(_grid.ind(i, j) - 1) = sources(_grid.ind(i, j) - 1) - epsilon * _assbl.Getlapdir(i, j) + epsilon * (_sigma(_grid.ind(i, j) - 1) - _sigma0(_grid.ind(i, j) - 1));
    //   }
    // }

    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      for (int i = 0; i < _grid.getNbPoints(0); i++)
      {

        for (int l = 0; l < ne; l++)
        {
          sources(_grid.ind(i, j) - 1) = sources(_grid.ind(i, j) - 1) - dsig(l, _grid.ind(i, j) - 1);
        }
        sources(_grid.ind(i, j) - 1) = sources(_grid.ind(i, j) - 1) - gradG(_grid.ind(i, j) - 1);
      }
    }

    _solver.compute(_assbl.GetMatrixdir());

    if (_solver.info() != Success)
    {
      std::cout << "" << std::endl;
      std::cout << "The decomposition not succeded" << std::endl;
    }
    _dir = _solver.solve(sources);
    if (_solver.info() != Success)
    {
      std::cout << "Solving the system suceeded, we have not a solution of the EIT problem!!" << std::endl;
    }

    _assbl.CLEAR();

    VectorXd indice(200);

    double deltasignorme = _dir.lpNorm<Eigen::Infinity>();
    _dir = _dir / deltasignorme;

    // Golden section serach.
    // t_moin = 0;
    // t_plus = 1;
    // a = (1 + sqrt(5.)) / 2;
    // t_l = t_moin + (t_plus - t_moin) / (1 + a);
    // t_r = t_plus + t_moin - t_l;
    // // stopping criterion.
    // //  while (abs(t_moin-t_plus) > 0.005){
    // while (abs(t_moin - t_plus) > 0.1)
    // {

    //   fl = F_TV(alpha, Imm, Umm, t_l, _sigma, _sigma0, _dir, angle_debut, angle_fin, Z);
    //   fr = F_TV(alpha, Imm, Umm, t_r, _sigma, _sigma0, _dir, angle_debut, angle_fin, Z);
    //   cout << "Fl:"
    //        << "  " << setprecision(17) << fl << "   " << t_l << endl;
    //   cout << "Fr:"
    //        << "  " << setprecision(17) << fr << "   " << t_r << endl;

    //   if (fl < fr)
    //   {
    //     t_plus = t_r;
    //     t_r = t_l;
    //     t_l = t_plus + t_moin - t_r;
    //   }
    //   else
    //   {
    //     t_moin = t_l;
    //     t_l = t_r;
    //     t_r = t_plus + t_moin - t_l;
    //   }
    // }

    // if (fl < fr)
    // {
    //   tt = t_r;
    //   Fsigma = fl;
    // }
    // else
    // {
    //   tt = t_l;
    //   Fsigma = fr;
    // }

    // std::cout<<deltasignorme<<std::endl;

    // for (double i = 0; i < 200 ; i++)
    // {
    //   indice(i) = (i)/100;
    //   cout << F(alpha,
    //             Imm,
    //             Umm,
    //             indice(i),
    //             _sigma,
    //             _sigma0,
    //             _dir,
    //             angle_debut,
    //             angle_fin,
    //             Z)
    //        << "  " << indice(i) << endl;
    // }

    // std::cout<<"Helooooooooooooooooooooooooooooooooooooooo"<<std::endl;

    // double deltasignorme = _dir.lpNorm<Eigen::Infinity>();
    // _dir = _dir / deltasignorme;

    // fl = F(alpha, Imm, Umm, 0.4, _sigma, _sigma0, _dir, angle_debut, angle_fin, Z);
    // Fsigma = fl;

    // for (int i = 0; i < _grid.getNbPoints(0); i++)
    // {
    //   for (int j = 0; j < _grid.getNbPoints(1); j++)
    //   {
    //     if (_grid.Levelset(i, j, alpha) <= 0)
    //     {
    //       _sigma(_grid.ind(i, j) - 1) = _sigma(_grid.ind(i, j) - 1) + tt * _dir(_grid.ind(i, j) - 1);
    //     }
    //   }
    // }

    //

   // double c = 1 * tt;

    // for (int l = 1; l < 5; l++)
    //   for (int i = 0; i < n; i++)
    //   {
    //     while (_sigma[i] < 0)
    //     {
    //       _sigma -= (l - 1) * 0.75 * tt * _dir;
    //       c = l * 0.75 * tt;
    //       _sigma += c * _dir;
    //     }
    //   }

    _sigma += 10 * _dir;

    if ((t % 5) == 0)
    {
      for (int i = 0; i < n; i++)
      {
        _sigma0[i] = _sigma[i];
      }
    }

    cout << setprecision(17) << "Fisigma ::"
         << "      " << Fsigma << "    " << tt << "     " << t << endl;

    _writer.write(_sigma, alpha, t);
  }
}

void Inverser::runcombinedinv(VectorXd &alpha,
                              VectorXd &angle_debut,
                              VectorXd &angle_fin)
{
  _grid.BuildT(alpha);
  _grid.inter(alpha);
  _grid.Zonelectrode(alpha, angle_debut, angle_fin);

  uint n = _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2);
  int ne = _grid.getNbPoints(2);

  std::cout << "Nx = " << _grid.getNbPoints(0)
            << ", Ny = " << _grid.getNbPoints(1)
            << std::endl;

  MatrixXd currents(ne, ne);
  MatrixXd measurments(ne, ne);
  MatrixXd U(ne, ne);
  MatrixXd u(ne, n);
  MatrixXd w(ne, n);
  VectorXd Im(ne), sources(n);
  MatrixXf dsig(ne, n);

  VectorXd sig(n),
      sig_terre(n);
  VectorXd contact(ne),
      contact_terre(ne);

  VectorXd descent_conductivity(n), descent_contact(ne);

  MatrixXd descent(ne, ne);
  double t_moin, t_plus, t_l, t_r, tt, a, fr, fl;
  double Fconductivity;

  // Initialisation des directions de descente
  descent_conductivity.setConstant(0);
  descent.setConstant(0);
  descent_contact.setConstant(0);

  // Initialisation de l'impedence et de la conductivite
  contact_terre.setConstant(1);
  // contact_terre(0) = 0.2;
  // contact_terre(5) = 2;
  sig_terre.setConstant(1);
  contact = contact_terre;
  sig = sig_terre;

  double epsilon = 0.001;

  // Ouvrir les fichier et faire la lecture
  std::string folder = _params->getString("output dir");
  for (int l = 0; l < ne - 1; l++)
  {
    std::string nomfichier = folder + "/measure-elec-" + std::to_string(l) + ".txt";
    std::ifstream fichier(nomfichier.c_str()); // le constructeur de std::ifstream n'accepte pas de string
    while (!fichier.eof())
    {
      for (int i = 0; i < ne; i++)
      {
        fichier >> setprecision(17) >> currents(l, i) >> measurments(l, i);
      }
    }

    for (int i = 0; i < ne; i++)
    {
      measurments(l, i) = generator(measurments(l, i));
    }
  }

  // The actual algorithm :

  for (int t = 0; t < 300; t++)
  {

    descent_contact.setConstant(0);
    descent_conductivity.setConstant(0);

    for (int i = 0; i < ne; i++)
    {
      descent_contact(i) = epsilon * (contact(i) - contact_terre(i));
    }

    for (uint l = 0; l < ne; l++)
    {

      for (int i = 0; i < ne; i++)
      {
        Im(i) = currents(l, i);
      }
      _pb->run(sig, alpha, Im, angle_debut, angle_fin, contact);
      for (int i = 0; i < n; i++)
      {
        u(l, i) = _pb->GetSolution()[i];
      }

      for (int i = 0; i < ne; i++)
      {
        U(l, i) = u(l, _grid.Get_nbsigne() + i);
        Im(i) = U(l, i) - measurments(l, i);
      }

      _pb->run(sig, alpha, Im, angle_debut, angle_fin, contact);
      for (int i = 0; i < n; i++)
      {
        w(l, i) = _pb->GetSolution()[i];
      }

      _assbl.BuildDirMatrix(alpha,
                            angle_debut,
                            angle_fin,
                            sig);
      _assbl.Buildsourcedir1(l,
                             u,
                             w,
                             alpha);
      _assbl.Buildsourcedir2(sig,
                             sig,
                             alpha);

      for (int i = 0; i < _grid.getNbPoints(0); i++)
      {
        for (int j = 0; j < _grid.getNbPoints(1); j++)
        {
          dsig(l, _grid.ind(i, j) - 1) = _assbl.Getdsigma(i, j);
        }
      }

      for (int i = 0; i < ne; i++)
      {
        descent(l, i) = _assbl.Compute_integral(u, w, alpha, angle_debut, angle_fin, i, l);
      }
    }

    // compute the decent direction for the impedence
    for (int i = 0; i < ne; i++)
    {
      for (int l = 0; l < ne; l++)
      {
        descent_contact(i) -= descent(l, i);
      }
    }

    // normalisation
    double norm = 0;
    for (int i = 0; i < ne; i++)
      norm += descent_contact(i) * descent_contact(i);
    norm = sqrt(norm);
    descent_contact = descent_contact / norm;
    descent_contact = 0.05 * descent_contact;

    sources.setConstant(0);
    for (int j = 0; j < _grid.getNbPoints(1); j++)
    {
      for (int i = 0; i < _grid.getNbPoints(0); i++)
      {
        for (int l = 0; l < ne; l++)
        {
          sources(_grid.ind(i, j) - 1) = sources(_grid.ind(i, j) - 1) - dsig(l, _grid.ind(i, j) - 1);
        }
        sources(_grid.ind(i, j) - 1) = sources(_grid.ind(i, j) - 1) - epsilon * _assbl.Getlapdir(i, j) + epsilon * (sig(_grid.ind(i, j) - 1) - sig_terre(_grid.ind(i, j) - 1));
      }
    }

    _solver.compute(_assbl.GetMatrixdir());
    if (_solver.info() != Success)
      std::cout << "The decomposition failed!" << std::endl;
    descent_conductivity = _solver.solve(sources);
    if (_solver.info() != Success)
      std::cout << "Solving the system  failed!" << std::endl;

    _assbl.CLEAR();

    norm = 0;
    for (int i = 0; i < n; i++)
      norm += descent_conductivity(i) * descent_conductivity(i);
    norm = sqrt(norm);
    descent_conductivity = descent_conductivity/norm ;

  

    //std::cout<<descent_conductivity<<std::endl;

    // Golden section serach.
    // t_moin = 0;
    // t_plus = 20;
    // a = (1 + sqrt(5.)) / 2;
    // t_l = t_moin + (t_plus - t_moin) / (1 + a);
    // t_r = t_plus + t_moin - t_l;
    // // stopping criterion.
    // // while (abs(t_moin-t_plus) > 0.005){
    // while (abs(t_moin - t_plus) > 0.05)
    // {

    //   fl = F(alpha,
    //          currents,
    //          measurments,
    //          t_l,
    //          sig,
    //          sig_terre,
    //          descent_conductivity,
    //          angle_debut,
    //          angle_fin,
    //          contact);
    //   fr = F(alpha,
    //          currents,
    //          measurments,
    //          t_r,
    //          sig,
    //          sig_terre,
    //          descent_conductivity,
    //          angle_debut,
    //          angle_fin,
    //          contact);
    //   cout << " Fl : "
    //        << "  " << fl << "   " << t_l << endl;
    //   cout << " Fr : "
    //        << "  " << fr << "   " << t_r << endl;

    //   if (fl < fr)
    //   {
    //     t_plus = t_r;
    //     t_r = t_l;
    //     t_l = t_plus + t_moin - t_r;
    //   }
    //   else
    //   {
    //     t_moin = t_l;
    //     t_l = t_r;
    //     t_r = t_plus + t_moin - t_l;
    //   }
    // }

    // if (fl < fr)
    // {
    //   tt = t_r;
    //   Fconductivity = fl;
    // }
    // else
    // {
    //   tt = t_l;
    //   Fconductivity = fr;
    // }

    for (int i = 0; i < _grid.getNbPoints(0); i++)
    {
      for (int j = 0; j < _grid.getNbPoints(1); j++)
      {
        if (_grid.Levelset(i, j, alpha) <= 0)
        {
          sig(_grid.ind(i, j) - 1) = sig(_grid.ind(i, j) - 1) + 10 * descent_conductivity(_grid.ind(i, j) - 1);
          //std::cout<<descent_conductivity(_grid.ind(i, j) - 1)<<std::endl;

        }
      }
    }

    for (int i = 0; i < ne; i++)
    {
      contact(i) = contact(i) - 10 * descent_contact(i);
      std::cout << contact(i)
                << "  " << descent_contact(i)
                << "  " << i << std::endl;
    }

    if ((t % 5) == 0)
    {
      sig_terre = sig;
      contact_terre = contact;
    }

    //Combined cost functions
    // double Fconductivity = F(alpha,
    //                          currents,
    //                          measurments,
    //                          100,
    //                          sig,
    //                          sig_terre,
    //                          descent_conductivity,
    //                          angle_debut,
    //                          angle_fin,
    //                          contact);

    //  VectorXd step(100);

    // for (uint i = 0; i <10 ; i++)
    // {
    //   step(i) = i / 1.;

    //   std::cout
    //       << cost_function_both(alpha,
    //                             angle_debut,
    //                             angle_fin,
    //                             sig,
    //                             sig_terre,
    //                             contact,
    //                             contact_terre,
    //                             currents,
    //                             measurments,
    //                             step(i),
    //                             descent_contact,
    //                             descent_conductivity)
    //   << "  "
    //   << step(i) << "  " << i << std::endl;
    // }

    // double Fconductivityandz = cost_function_both(alpha,
    //                                               angle_debut,
    //                                               angle_fin,
    //                                               sig,
    //                                               sig_terre,
    //                                               contact,
    //                                               contact_terre,
    //                                               currents,
    //                                               measurments,
    //                                               0.05,
    //                                               descent_contact,
    //                                               descent_conductivity);
    double Fconductivity = F(alpha,
                             currents,
                             measurments,
                             10,
                             sig,
                             sig_terre,
                             descent_conductivity,
                             angle_debut,
                             angle_fin,
                             contact);

    double Fimpedance = cost_function(alpha,
                                angle_debut,
                                angle_fin,
                                sig,
                                contact,
                                contact_terre,
                                currents,
                                measurments,
                                10,
                                descent_contact);

    cout << setprecision(17) << "Fisigma ::"
         << "      " << Fconductivity << "     " << t << endl;
    cout << setprecision(17) << "Fcontact ::"
         << "      " << Fimpedance << "     " << t << endl;
    // cout << setprecision(17) << "Fcontact ::"
    //      << "      " <<  F(alpha,
    //                          currents,
    //                          measurments,
    //                          10,
    //                          sig,
    //                          sig_terre,
    //                          descent_conductivity,
    //                          angle_debut,
    //                          angle_fin,
    //                          contact) << "     " << t << endl;

    _writer.write(sig, alpha, t);
  }
}