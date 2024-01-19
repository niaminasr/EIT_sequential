#include "heatProblem.hpp"
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

HeatProblem::HeatProblem(Parameters &p)
{

  _params = &p;
  _k = p.getReal("diffusion coefficient", 1.);
  _grid = Grid(p);
  _assbl = EITAssembler(&_grid, this);
  _writer = Writer(p, _grid);

  squaresize = 2.;

  _x.resize(_grid.getNbPoints(0));
  for (int i = 0; i < _grid.getNbPoints(0); ++i)
    _x[i] = ((i + 1) - 0.5) * _grid.getSpacing(0) - squaresize;

  _y.resize(_grid.getNbPoints(1));
  for (int i = 0; i < _grid.getNbPoints(1); ++i)
    _y[i] = ((i + 1) - 0.5) * _grid.getSpacing(1) - squaresize;
}

void HeatProblem::run(VectorXd _sigma,
                      VectorXd &alpha,
                      VectorXd &Imm,
                      VectorXd &angle_debut,
                      VectorXd &angle_fin,
                      VectorXd &Z)
{

  _assbl.BuildLaplacianMatrix(alpha, angle_debut, angle_fin, _sigma,Z);
  _assbl.BuildSourceTerm(alpha, Imm, Z);

  _sol.resize(_grid.getNbPoints(0) * _grid.getNbPoints(1) +
              _grid.Get_nbsigne() + _grid.getNbPoints(2));

  _sol.setConstant(0.0);

  auto &solver = _solver;
  auto A = _assbl.GetMatrix();
  auto b = _assbl.GetSourceTerm();
  solver.compute(A);
  _sol = solver.solve(b);

  double sumu = 0;
  for (int i = 0; i < _grid.getNbPoints(2); i++)
  {
    sumu = sumu + _sol[_grid.Get_nbsigne() + i];
  }

  for (int i = 0; i < _grid.getNbPoints(0) * _grid.getNbPoints(1) +
                          _grid.Get_nbsigne() + _grid.getNbPoints(2);
       i++)
  {
    _sol[i] -= sumu / _grid.getNbPoints(2);
  }
  _writer.write(_sol, alpha, 997);

  _assbl.CLEAR();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HeatProblem::runwithsource(Eigen::VectorXd &_sigma,
                                Eigen::VectorXd &alpha,
                                Eigen::VectorXd &angle_debut,
                                Eigen::VectorXd &angle_fin,
                                VectorXd &Z)
{

  float temps;
  clock_t t1, t2;

  _grid.BuildT(alpha);
  _grid.DLevelset(alpha);
  _grid.inter(alpha);
  _assbl.BuildLaplacianMatrix(alpha, angle_debut, angle_fin, _sigma, Z);
  _assbl.BuildSourceTermMAN(alpha, Z);

  _sol.resize(_grid.getNbPoints(0) * _grid.getNbPoints(1) +
              _grid.Get_nbsigne() + _grid.getNbPoints(2));
  yy.resize(_grid.getNbPoints(0) * _grid.getNbPoints(1) +
            _grid.Get_nbsigne() + _grid.getNbPoints(2));
  xx.resize(_grid.getNbPoints(0) * _grid.getNbPoints(1) +
            _grid.Get_nbsigne() + _grid.getNbPoints(2));
  residu.resize(_grid.getNbPoints(0) * _grid.getNbPoints(1) +
                _grid.Get_nbsigne() + _grid.getNbPoints(2));

  _sol.setConstant(0.0);

  //////////////////////////////////////////
  float time = 0;
  t1 = clock();
  _solver.compute(_assbl.GetMatrix());
  if (_solver.info() != Success)
  {
    std::cout << "the decomposition of the matrix failed" << std::endl;
    return;
  }
  if (_solver.info() == Success)
  {
    std::cout << "the decomposition succeded" << std::endl;
  }
  _sol = _solver.solve(_assbl.GetSourceTermMan());
  if (_solver.info() == Success)
  {
    std::cout << "solving the system succeded" << std::endl;
  }
  t2 = clock();
  temps = (float)(t2 - t1) / CLOCKS_PER_SEC;
  printf("temps = %f\n", temps);
  /////////////////////////////////////////
  double pond1, pond2;
  pond1 = _sol[_grid.ind(_grid.getNbPoints(0) / 2 - 1,
                         _grid.getNbPoints(1) / 2 - 1)];

  pond2 = _grid.cond_inter(_x[_grid.getNbPoints(0) / 2 - 1],
                           _y[_grid.getNbPoints(1) / 2 - 1]);

  for (int i = 0; i < _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2); i++)
  {
    _sol[i] = _sol[i] - pond1 + pond2;
  }
  ///////////////////////////////////////

  _dxsol.resize(_grid.getNbPoints(0), _grid.getNbPoints(1));

  for (unsigned short i = 0; i < ((unsigned short)(_grid.getNbPoints(0))); i++)
  {
    for (unsigned short j = 0; j < ((unsigned short)(_grid.getNbPoints(1))); j++)
    {

      if (_grid.Levelset(i, j, alpha) <= 0)
      {

        uup = _sol[_grid.ind(i + 1, j) - 1];
        xp = _x[i + 1];
        if ((_grid.get_tx2(i, j)) != 0)
        {
          uup = _sol[_grid.get_tx2(i, j) - 1];
          xp = _grid.get_xinter(i, j);
        }

        uum = _sol[_grid.ind(i - 1, j) - 1];
        xm = _x[i - 1];
        if ((_grid.get_tx2(i - 1, j)) != 0)
        {
          uum = _sol[_grid.get_tx2(i - 1, j) - 1];
          xm = _grid.get_xinter(i - 1, j);
        }

        _dxsol(i, j) = (uup - uum) / (xp - xm);
      }
    }
  }

  /////////////////////////////////////////
  for (int i = 0; i < _grid.getNbPoints(0) * _grid.getNbPoints(1) + _grid.Get_nbsigne() + _grid.getNbPoints(2); i++)
  {
    residu[i] = 0;
  }
  residu = _assbl.GetMatrix() * _sol - _assbl.GetSourceTermMan();
  errmax = 0;
  errdxmax = 0.;
  for (unsigned short j = 0; j < ((unsigned short)(_grid.getNbPoints(1))); j++)
  {
    for (unsigned short i = 0; i < ((unsigned short)(_grid.getNbPoints(0))); i++)
    {

      exactsol = 0.;
      if (_grid.Levelset(i, j, alpha) <= 0)
        exactsol = _grid.cond_inter(_x[i], _y[j]);
      if (_grid.Levelset(i, j, alpha) <= 0)
      {
        dxsolexact = _sigma(_grid.ind(i,j))*_y[j] * cos(_x[i] * _y[j]);
        if (errmax < abs((_sol(_grid.ind(i, j) - 1) - exactsol)))
        {
          errmax = abs((_sol(_grid.ind(i, j) - 1) - exactsol));
        }

        if (errdxmax < abs(_dxsol(i, j) - dxsolexact))
        {
          errdxmax = abs(_dxsol(i, j) - dxsolexact);
        }
      }
    }
  }

  resmax = 0.;
  for (unsigned short i = 0; i < ((unsigned short)(_grid.getNbPoints(0))); i++)
  {
    for (unsigned short j = 0; j < ((unsigned short)(_grid.getNbPoints(1))); j++)
    {
      if (_grid.Levelset(i, j, alpha) < 0.)
      {
        if (resmax < abs(residu(_grid.ind(i, j) - 1)))
        {
          resmax = abs(residu(_grid.ind(i, j) - 1));
        }
      }
    }
  }

  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "errmax:"
            << " " << std::setprecision(17) << _grid.getNbPoints(0) << "   "
            << errmax << "   " << std::setprecision(20)
            << _grid.getNbPoints(0) * errmax << std::endl;
  std::cout << "errdxmax:"
            << " " << std::setprecision(17) << _grid.getNbPoints(0) << "   "
            << errdxmax << "   " << std::setprecision(20)
            << _grid.getNbPoints(0) * errdxmax << std::endl;
  std::cout << "résidus:"
            << " " << std::setprecision(20) << resmax << std::endl;

  t1 = clock();
  _writer.write(_sol, alpha, 0);
  _assbl.CLEAR();
  t2 = clock();
  temps = (float)(t2 - t1) / CLOCKS_PER_SEC;
  printf("temps = %f\n", temps);
}
