#include "inverseforme.hpp"
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <complex>
#include "Eigen/Eigen/Eigen"
#include "Eigen/Eigen/Dense"
#include "Eigen/Eigen/Core"
#include "Eigen/Eigen/Sparse"
#include "Eigen/Eigen/IterativeLinearSolvers"
#include "Eigen/unsupported/Eigen/src/IterativeSolvers/GMRES.h"

using namespace std;

Inverserr::Inverserr(Parameters &p, HeatProblem &pb) : _pb(&pb), _params(&p)
{

  _grid = Grid(p);
  _assbl = EITAssembler(&_grid, _pb);        // not sure about _pb
  _assblf = HeatAssemblerforme(&_grid, _pb); // not sure about _pb
  _writer = Writer(p, _grid);

  // This also calls copy constructor

  squaresize = 2.;

  _x.resize(_grid.getNbPoints(0)); // (xmin+h, ..., xmax-h)
  for (int i = 0; i < _grid.getNbPoints(0); ++i)
    _x[i] = ((i + 1) - 0.5) * _grid.getSpacing(0) - squaresize;

  _y.resize(_grid.getNbPoints(1));
  for (int i = 0; i < _grid.getNbPoints(1); ++i)
    _y[i] = ((i + 1) - 0.5) * _grid.getSpacing(1) - squaresize;

  dx = _grid.getSpacing(0);
  dy = _grid.getSpacing(1);
}

void Inverserr::run_inv_cond_in_any_geom()
{
  int i, j, l, k, ind_g, m;
  double d, verif, F;
  double epsilon = 0.001;
  double pi = acos(-1.);
  int ne = _grid.getNbPoints(2);

  Inverser _inver(*_params, *_pb);
  Eigen::VectorXd sigma, Z(ne), Z_terre(ne);
  Eigen::VectorXd Im, norme, norme_n;
  Eigen::VectorXd _alpha, _alpha_0, Theta_begin, Theta_begin_0, Theta_end, Theta_end_0;
  _alpha_0.resize(7);
  _alpha.resize(7);
  Theta_begin.resize(ne);
  Theta_begin_0.resize(ne);
  Theta_end.resize(ne);
  Theta_end_0.resize(ne);

  Im.resize(ne);

  MatrixXf Imm(ne, ne);
  MatrixXf Umm(ne, ne);
  MatrixXf aa(ne, ne);
  MatrixXf bb(ne, ne);

  _alpha.setConstant(0);

  // int shape = _params->getInt("shape");
  // if (shape == 1)
  // {
  //   _alpha[0] = 1.4;
  //   _alpha[1] = 0.0;
  //   _alpha[2] = 0.41;
  //   _alpha[3] = 0.0;
  //   _alpha[4] = 0.1;
  //   _alpha[5] = 0.0;
  //   _alpha[6] = 0.;
  // }
  // else if (shape == 2)
  // {
  //   _alpha[0] = 1.51;
  //   _alpha[1] = 0.01;
  //   _alpha[2] = 0.05;
  //   _alpha[3] = 0.2;
  //   _alpha[4] = 0.035;
  //   _alpha[5] = 0.01;
  //   _alpha[6] = 0.1;
  // }
  // else if (shape == 3)
  // {
  //   _alpha[0] = 1.6;
  //   _alpha[1] = 0.002;
  //   _alpha[2] = 0.01;
  //   _alpha[3] = 0.003;
  //   _alpha[4] = 0.035;
  //   _alpha[5] = 0.2;
  //   _alpha[6] = 0.15;
  // }
  int shape = _params->getInt("shape");
  if (shape == 1)
  {
    _alpha[0] = 1.5;
    _alpha[1] = 0.0;
    _alpha[2] = 0.0;
    _alpha[3] = 0.0;
    _alpha[4] = 0.0;
    _alpha[5] = 0.0;
    _alpha[6] = 0.0;
  }
  else if (shape == 2)
  {
    _alpha[0] = 1.51;
    _alpha[1] = 0.01;
    _alpha[2] = 0.05;
    _alpha[3] = 0.2;
    _alpha[4] = 0.035;
    _alpha[5] = 0.01;
    _alpha[6] = 0.1;
  }
  else if (shape == 3)
  {
    _alpha[0] = 1.6;
    _alpha[1] = 0.002;
    _alpha[2] = 0.01;
    _alpha[3] = 0.003;
    _alpha[4] = 0.035;
    _alpha[5] = 0.2;
    _alpha[6] = 0.15;
  }

  for (i = 0; i < Theta_begin.size(); i++)
  {
    Theta_begin(i) = ((i) * 2. * pi / ne - pi);
  }

  double theta_plus = _params->getReal("theta_plus");
  for (i = 0; i < Theta_end.size(); i++)
  {
    Theta_end(i) = _grid.angleplus(0.35, Theta_begin(i), _alpha);
  }

  for (i = 0; i < Theta_begin.size(); i++)
  {
    cout << setprecision(17) << Theta_begin(i) << "  " << Theta_end(i) << endl;
  }

  _grid.BuildT(_alpha);
  _grid.inter(_alpha);
  _grid.DLevelset(_alpha);
  _grid.Zonelectrode(_alpha, Theta_begin, Theta_end);
  _grid.ordre();
  _grid.coord_int();
  uint n = _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2);
  MatrixXd u;
  u.resize(ne, n);
  sigma.resize(n);

  sigma.setConstant(1);
  Z.setConstant(1);

  int testid = _params->getInt("testid");
  //   if ((pow(_x[i]+0.6,2)+pow(_y[j]-0.4,2)) <= 0.250)
  //         sigma(_grid.ind(i,j))=0.2;

  // ***Test_2***:: add inclusion in the center
  // if ((pow(_x[i],2)+pow(_y[j],2)) <= 0.3)
  //         sigma(_grid.ind(i,j))=2;

  // ***Test_3***:: add two inclusions
  // if ((pow(_x[i]-0.2,2)+pow(_y[j]+0.5,2)) <= 0.3)
  //         sigma(_grid.ind(i,j))=1.5;
  // if ((pow(_x[i]+0.6,2)+pow(_y[j]-0.7,2)) <= 0.10
  //         sigma(_grid.ind(i,j))=0.5;

  if (shape == 1)
  {
    if (testid == 1)
    {
      // VectorXd alphas_2(7), alphas_3(7);
      // alphas_2.setConstant(0);
      // alphas_2[0] = 0.6;
      // alphas_2[1] = -0.005;
      // alphas_2[2] = -0.025;
      // alphas_2[3] = -0.0075;
      // alphas_2[4] = -0.0007;
      // alphas_2[5] = -0.09;
      // alphas_2[6] = -0.03;

      // alphas_3.setConstant(0);
      // alphas_3[0] = 0.6;
      // alphas_3[1] = 0.005;
      // alphas_3[2] = 0.025;
      // alphas_3[3] = 0.0075;
      // alphas_3[4] = 0.0007;
      // alphas_3[5] = 0.09;
      // alphas_3[6] = 0.03;
      //  for (i = 0; i < _grid.getNbPoints(0); i++)
      //   for (j = 0; j < _grid.getNbPoints(1); j++)
      //   {
      //     if (_grid.f(_x[i], _y[j], -0.9, 0.081, alphas_2) <= 0)
      //       sigma(_grid.ind(i, j)) = 0.42;

      //     if (_grid.f(_x[i], _y[j], 0.9, 0.0, alphas_3) <= 0)
      //       sigma(_grid.ind(i, j)) = 0.42;

      //   }

      Z.setConstant(2);

      for (i = 0; i < _grid.getNbPoints(0); i++)
        for (j = 0; j < _grid.getNbPoints(1); j++)
        {
          if ((pow(_x[i] + 0.6, 2) + pow(_y[j] - 0.4, 2)) <= 0.250)
            sigma(_grid.ind(i, j)) = 0.2;
        }

      _writer.write(sigma, _alpha, 999);
    }
    else if (testid == 2)
    {

      for (i = 0; i < _grid.getNbPoints(0); i++)
        for (j = 0; j < _grid.getNbPoints(1); j++)
        {
          if ((pow(_x[i], 2) + pow(_y[j], 2)) <= 0.3)
            sigma(_grid.ind(i, j)) = 2;
        }
       Z.setConstant(0.8);
      // VectorXd alphas_1(7), alphas_2(7), alphas_3(7);
      // alphas_1.setConstant(0);
      // alphas_1[0] = 0.3;
      // alphas_1[3] = 0.018;
      // alphas_1[5] = 0.056;
      // alphas_1[6] = 0.018;

      // alphas_2.setConstant(0);
      // alphas_2[0] = 0.6;
      // alphas_2[1] = -0.005;
      // alphas_2[2] = -0.025;
      // alphas_2[3] = -0.0075;
      // alphas_2[4] = -0.0007;
      // alphas_2[5] = -0.09;
      // alphas_2[6] = -0.03;

      // alphas_3.setConstant(0);
      // alphas_3[0] = 0.6;
      // alphas_3[1] = 0.005;
      // alphas_3[2] = 0.025;
      // alphas_3[3] = 0.0075;
      // alphas_3[4] = 0.0007;
      // alphas_3[5] = 0.09;
      // alphas_3[6] = 0.03;
      // for (i = 0; i < _grid.getNbPoints(0); i++)
      //   for (j = 0; j < _grid.getNbPoints(1); j++)
      //   {
      //     if (_grid.f(_x[i], _y[j], 0.1, 0.7, alphas_1) <= 0)
      //       sigma(_grid.ind(i, j)) = 3;

      //     if (_grid.f(_x[i], _y[j], -1, 0.081, alphas_2) <= 0)
      //       sigma(_grid.ind(i, j)) = 0.42;

      //     if (_grid.f(_x[i], _y[j], 0.9, 0.0, alphas_3) <= 0)
      //       sigma(_grid.ind(i, j)) = 0.42;
      //   }
      _writer.write(sigma, _alpha, 999);
    }
    else if (testid == 3)
    {
      Z.setConstant(1);
      for (i = 0; i < _grid.getNbPoints(0); i++)
      {
        for (j = 0; j < _grid.getNbPoints(1); j++)
        {
          if ((pow(_x[i] - 0.8, 2) + pow(_y[j] - 0.3, 2)) <= 0.2)
            sigma(_grid.ind(i, j)) = 1.5;
          if ((pow(_x[i] + 0.3, 2) + pow(_y[j] + 0.8, 2)) <= 0.15)
            sigma(_grid.ind(i, j)) = 0.2;
        }
      }
      // VectorXd alphas_1(7), alphas_2(7), alphas_3(7), alphas_4(7);
      // alphas_1.setConstant(0);
      // alphas_1[0] = 0.3;
      // alphas_1[3] = 0.018;
      // alphas_1[5] = 0.056;
      // alphas_1[6] = 0.018;

      //   alphas_2.setConstant(0);
      // alphas_2[0] = 0.6;
      // alphas_2[1] = -0.005;
      // alphas_2[2] = -0.025;
      // alphas_2[3] = -0.0075;
      // alphas_2[4] = -0.0007;
      // alphas_2[5] = -0.09;
      // alphas_2[6] = -0.03;

      // alphas_3.setConstant(0);
      // alphas_3[0] = 0.6;
      // alphas_3[1] = 0.005;
      // alphas_3[2] = 0.025;
      // alphas_3[3] = 0.0075;
      // alphas_3[4] = 0.0007;
      // alphas_3[5] = 0.09;
      // alphas_3[6] = 0.03;

      // alphas_4.setConstant(0);
      // alphas_4[0] = 0.2;
      // alphas_4[6] = 0.03;
      // for (i = 0; i < _grid.getNbPoints(0); i++)
      // {
      //   for (j = 0; j < _grid.getNbPoints(1); j++)
      //   {

      //   if (_grid.f(_x[i], _y[j], 0.1, 0.7, alphas_1) <= 0)
      //       sigma(_grid.ind(i, j)) = 2.5;

      //     if (_grid.f(_x[i], _y[j], -0.9, 0.081, alphas_2) <= 0)
      //       sigma(_grid.ind(i, j)) = 0.42;

      //     if (_grid.f(_x[i], _y[j], 0.9, 0.0, alphas_3) <= 0)
      //       sigma(_grid.ind(i, j)) = 0.42;

      //     if (_grid.f(_x[i], _y[j], 0.0, -0.5, alphas_4) <= 0)
      //       sigma(_grid.ind(i, j)) = 0.06;
      //   }
      // }
      _writer.write(sigma, _alpha, 999);
    }
  }
  else if (shape == 2)
  {
    if (testid == 1)
    {
      Z.setConstant(2);
      Z(6) = 2.5;
      for (i = 0; i < _grid.getNbPoints(0); i++)
        for (j = 0; j < _grid.getNbPoints(1); j++)
        {
          if ((pow(_x[i] + 0.3, 2) + pow(_y[j] + 0.8, 2)) <= 0.3)
            sigma(_grid.ind(i, j)) = 0.2;
        }
      _writer.write(sigma, _alpha, 999);
    }
    else if (testid == 2)
    {
      Z.setConstant(3);
      for (i = 0; i < _grid.getNbPoints(0); i++)
        for (j = 0; j < _grid.getNbPoints(1); j++)
        {
          if ((pow(_x[i], 2) + pow(_y[j], 2)) <= 0.25)
            sigma(_grid.ind(i, j)) = 0.2;
        }

      _writer.write(sigma, _alpha, 999);
    }
    else if (testid == 3)
    {
      Z.setConstant(1.5);
      Z(6) = 3;
      Z(8) = 3;
      for (i = 0; i < _grid.getNbPoints(0); i++)
      {
        for (j = 0; j < _grid.getNbPoints(1); j++)
        {
          if ((pow(_x[i] - 0.8, 2) + pow(_y[j] - 0.3, 2)) <= 0.2)
            sigma(_grid.ind(i, j)) = 2;
          if ((pow(_x[i] + 0.3, 2) + pow(_y[j] + 0.8, 2)) <= 0.15)
            sigma(_grid.ind(i, j)) = 0.2;
        }
      }
      _writer.write(sigma, _alpha, 999);
    }
  }
  else if (shape == 3)
  {
    if (testid == 1)
    {
      Z.setConstant(0.7);
      for (i = 0; i < _grid.getNbPoints(0); i++)
        for (j = 0; j < _grid.getNbPoints(1); j++)
        {
          if ((pow(_x[i] - 0.9, 2) + pow(_y[j] - 0.3, 2)) <= 0.25)
            sigma(_grid.ind(i, j)) = 0.1;
        }
      _writer.write(sigma, _alpha, 999);
    }
    else if (testid == 2)
    {
      Z.setConstant(1.5);
      Z(2) = 0.8;
      for (i = 0; i < _grid.getNbPoints(0); i++)
        for (j = 0; j < _grid.getNbPoints(1); j++)
        {
          if ((pow(_x[i], 2) + pow(_y[j], 2)) <= 0.3)
            sigma(_grid.ind(i, j)) = 2;
        }
      _writer.write(sigma, _alpha, 999);
    }
    else if (testid == 3)
    {
      Z.setConstant(2);
      Z(1) = 0.8;
      Z(12) = 0.9;
      for (i = 0; i < _grid.getNbPoints(0); i++)
      {
        for (j = 0; j < _grid.getNbPoints(1); j++)
        {
          if ((pow(_x[i] - 1, 2) + pow(_y[j] - 0.5, 2)) <= 0.25)
            sigma(_grid.ind(i, j)) = 1.7;
          if ((pow(_x[i] + 0.4, 2) + pow(_y[j] + 0.6, 2)) <= 0.2)
            sigma(_grid.ind(i, j)) = 0.2;
        }
      }
    }
    _writer.write(sigma, _alpha, 999);
  }
  Imm.diagonal<0>().setConstant(1.4142);
  Imm.diagonal<1>().setConstant(-1.4142);

  std::string folder = _params->getString("output dir");
  for (int l = 0; l < ne - 1; l++)
  {

    for (int i = 0; i < ne; i++)
    {
      Im(i) = Imm(l, i);
    }

    _pb->run(sigma, _alpha, Im, Theta_begin, Theta_end, Z);
    u.row(l) = _pb->GetSolution();

    cout << "\nwrite measure file for test " << l << endl;

    std::string nomfichier = folder + "/measure-elec-" + std::to_string(l) + ".txt";
    std::ofstream fichier(nomfichier.c_str()); // le constructeur de std::ofstream n'accepte pas de string
    for (int i = 0; i < ne; i++)
    {
      cout << u(l, _grid.Get_nbsigne() + i) << "  " << Im(i) << endl;
      fichier << std::setprecision(7) << Im(i) << " " << u(l, _grid.Get_nbsigne() + i) << std::endl;
    }
    fichier.close();
    cout << "done for test " << l << endl;
  }

  Z.setConstant(1);
  sigma.setConstant(1);

  //_inver.runinvcontact(_alpha, Theta_begin, Theta_end, Z);
  _inver.runcombinedinv(_alpha, Theta_begin, Theta_end);
}