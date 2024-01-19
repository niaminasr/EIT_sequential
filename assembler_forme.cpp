#include "assembler_forme.hpp"
#include "heatProblem.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "Eigen/Eigen/Eigen"
#include"Eigen/Eigen/Dense"
#include "Eigen/Eigen/Sparse"

using namespace std;
using namespace Eigen;

HeatAssemblerforme::HeatAssemblerforme(Grid* g, HeatProblem* pb) {

  _grid = g;
  _pb   = pb;
  mat1.resize(8*_grid->getNbPoints(0)*_grid->getNbPoints(1) + 10* _grid->Get_nbsigne());
  mat2.resize(8*_grid->getNbPoints(0)*_grid->getNbPoints(1) + 10* _grid->Get_nbsigne());
  mat3.resize(8*_grid->getNbPoints(0)*_grid->getNbPoints(1) + 10* _grid->Get_nbsigne());
  Im.resize(_grid->getNbPoints(2));
  Um.resize(_grid->getNbPoints(2));
  Nx = _grid->getNbPoints(0);
  Ny = _grid->getNbPoints(1);
  Ne = _grid->getNbPoints(2);
  dx = _grid->getSpacing(0);
  dy = _grid->getSpacing(1);
  
  squaresize = 2.;

  _x.resize(Nx); // (xmin+h, ..., xmax-h)
  for (int i = 0 ; i < Nx ; ++i)
  _x[i] =  ((i+1)-0.5)*dx - squaresize;

  _y.resize(Ny);
  for (int i = 0 ; i < Ny ; ++i)
  _y[i] =  ((i+1)-0.5)*dy - squaresize;
  _eta = 0.0005;
}

void HeatAssemblerforme::fluxx_forme(int i, int j){

  pi = acos(-1.);

  fix=0;
  fiy=0;

  //calcul normale sur le point d'interface
  fix = (_grid->get_u(i,j)*(_x[i+1]-_grid->get_xinter(i,j))-_grid->get_u(i+1,j)*(_x[i]-_grid->get_xinter(i,j)))/dx;
  fiy = (_grid->get_v(i,j)*(_x[i+1]-_grid->get_xinter(i,j))-_grid->get_v(i+1,j)*(_x[i]-_grid->get_xinter(i,j)))/dx;
  fiix = fix/sqrt(pow(fix,2)+pow(fiy,2));
  fiiy = fiy/sqrt(pow(fix,2)+pow(fiy,2));
  fix = fiix;
  fiy = fiiy;

  
  // cout<<M << "\n"<<endl;
 
    Eigen:: MatrixXd M, M_ind;
    M.resize(6,2);
    M_ind.resize(6,2);


    M(0,0) = _x[i];
    M(0,1) = _y[j+1];
    M(1,0) = _x[i+1];
    M(1,1) = _y[j+1];
    M(2,0) = _x[i];
    M(2,1) = _y[j];
    M(3,0) = _x[i+1];
    M(3,1) = _y[j];
    M(4,0) = _x[i];
    M(4,1) = _y[j-1];
    M(5,0) = _x[i+1];
    M(5,1) = _y[j-1];

    M_ind(0,0) = i;
    M_ind(0,1) = j+1;
    M_ind(1,0) = i+1;
    M_ind(1,1) = j+1;
    M_ind(2,0) = i;
    M_ind(2,1) = j;
    M_ind(3,0) = i+1;
    M_ind(3,1) = j;
    M_ind(4,0) = i;
    M_ind(4,1) = j-1;
    M_ind(5,0) = i+1;
    M_ind(5,1) = j-1;

    xj = _grid->get_xinter(i,j);
    yj = _y[j];
    _grid->discret_scheme(i,j,fix,fiy,xj,yj,M, M_ind);

    xk =_grid->get_xk(i,j);
    yk =_grid->get_yk(i,j);
    indk = _grid->get_indk(i,j);

    xi =_grid->get_xi(i,j);
    yi =_grid->get_yi(i,j);
    indi = _grid->get_indi(i,j);

    // calcul des coef du gradient (interpolation lineaire sur les 3 points)
    denom = (xj-xk)*(yj-yi)-(xj-xi)*(yj-yk);
    alphaj = (yk-yi)/denom;
    betaj = (xi-xk)/denom;

    denom = (xi-xk)*(yi-yj)-(xi-xj)*(yi-yk);
    alphai = (yk-yj)/denom;
    betai = (xj-xk)/denom;

    denom = (xk-xj)*(yk-yi)-(xk-xi)*(yk-yj);
    alphak = (yj-yi)/denom;
    betak = (xi-xj)/denom;

    // indk
    mat1[m] = _grid->get_tx2(i,j)-1;
    mat2[m] = indk-1;
    mat3[m] = _grid->coef(_grid->get_xinter(i,j),_y[j])*(alphak*fix + betak*fiy);
    
    m = m +1;

    // indi
    mat1[m] = _grid->get_tx2(i,j)-1;
    mat2[m] = indi-1;
    mat3[m] = _grid->coef(_grid->get_xinter(i,j),_y[j])*(alphai*fix + betai*fiy);
    m = m +1;

    // int(i,i+1),j
    mat1[m] = _grid->get_tx2(i,j)-1;
    mat2[m] = _grid->get_tx2(i,j)-1;
    mat3[m] = _grid->coef(_grid->get_xinter(i,j),_y[j])*(alphaj*fix + betaj*fiy);
    m = m +1;
}


void HeatAssemblerforme::fluxy_forme(int i, int j){

  pi = acos(-1.);

  fix=0;
  fiy=0;

  //calcul normale sur le point d'interface
  fix = (_grid->get_u(i,j)*(_y[j+1]-_grid->get_yinter(i,j))-_grid->get_u(i,j+1)*(_y[j]-_grid->get_yinter(i,j)))/dy;
  fiy = (_grid->get_v(i,j)*(_y[j+1]-_grid->get_yinter(i,j))-_grid->get_v(i,j+1)*(_y[j]-_grid->get_yinter(i,j)))/dy;
  fiix = fix/sqrt(pow(fix,2)+pow(fiy,2));
  fiiy = fiy/sqrt(pow(fix,2)+pow(fiy,2));
  fix = fiix;
  fiy = fiiy;


    Eigen:: MatrixXd M, M_ind;
    M.resize(6,2);
    M_ind.resize(6,2);

    M(0,0) = _x[i-1];
    M(0,1) = _y[j+1];
    M(1,0) = _x[i-1];
    M(1,1) = _y[j];
    M(2,0) = _x[i];
    M(2,1) = _y[j+1];
    M(3,0) = _x[i];
    M(3,1) = _y[j];
    M(4,0) = _x[i+1];
    M(4,1) = _y[j+1];
    M(5,0) = _x[i+1];
    M(5,1) = _y[j];

    M_ind(0,0) = i-1;
    M_ind(0,1) = j+1;
    M_ind(1,0) = i-1;
    M_ind(1,1) = j;
    M_ind(2,0) = i;
    M_ind(2,1) = j+1;
    M_ind(3,0) = i;
    M_ind(3,1) = j;
    M_ind(4,0) = i+1;
    M_ind(4,1) = j+1;
    M_ind(5,0) = i+1;
    M_ind(5,1) = j;

    xj = _x[i];
    yj = _grid->get_yinter(i,j);
    _grid->discret_scheme(i,j,fix,fiy,xj,yj,M, M_ind);

    xk =_grid->get_xk(i,j);
    yk =_grid->get_yk(i,j);
    indk = _grid->get_indk(i,j);

    xi =_grid->get_xi(i,j);
    yi =_grid->get_yi(i,j);
    indi = _grid->get_indi(i,j);
    
    // calcul des coef du gradient (interpolation lineaire sur les 3 points)
    denom = (xj-xk)*(yj-yi)-(xj-xi)*(yj-yk);
    alphaj = (yk-yi)/denom;
    betaj = (xi-xk)/denom;

    denom = (xi-xk)*(yi-yj)-(xi-xj)*(yi-yk);
    alphai = (yk-yj)/denom;
    betai = (xj-xk)/denom;

    denom = (xk-xj)*(yk-yi)-(xk-xi)*(yk-yj);
    alphak = (yj-yi)/denom;
    betak = (xi-xj)/denom;

    // indk
    mat1[m] = _grid->get_ty2(i,j)-1;
    mat2[m] = indk-1;
    mat3[m] = _grid->coef(_x[i],_grid->get_yinter(i,j))*(alphak*fix + betak*fiy);
    m = m +1;

    // indi
    mat1[m] = _grid->get_ty2(i,j)-1;
    mat2[m] = indi-1;
    mat3[m] = _grid->coef(_x[i],_grid->get_yinter(i,j))*(alphai*fix + betai*fiy);
    m = m +1;

    // int(i,i+1),j
    mat1[m] = _grid->get_ty2(i,j)-1;
    mat2[m] = _grid->get_ty2(i,j)-1;
    mat3[m] = _grid->coef(_x[i],_grid->get_yinter(i,j))*(alphaj*fix + betaj*fiy);
    m = m +1;

}



void HeatAssemblerforme::electrodefluxx_forme(int i, int j){

  pi = acos(-1.);

  fix=0;
  fiy=0;

  //calcul normale sur le point d'interface
  fix = (_grid->get_u(i,j)*(_x[i+1]-_grid->get_xinter(i,j))-_grid->get_u(i+1,j)*(_x[i]-_grid->get_xinter(i,j)))/dx;
  fiy = (_grid->get_v(i,j)*(_x[i+1]-_grid->get_xinter(i,j))-_grid->get_v(i+1,j)*(_x[i]-_grid->get_xinter(i,j)))/dx;
  fiix = fix/sqrt(pow(fix,2)+pow(fiy,2));
  fiiy = fiy/sqrt(pow(fix,2)+pow(fiy,2));
  fix = fiix;
  fiy = fiiy;

    Eigen:: MatrixXd M, M_ind;
    M.resize(6,2);
    M_ind.resize(6,2);


    M(0,0) = _x[i];
    M(0,1) = _y[j+1];
    M(1,0) = _x[i+1];
    M(1,1) = _y[j+1];
    M(2,0) = _x[i];
    M(2,1) = _y[j];
    M(3,0) = _x[i+1];
    M(3,1) = _y[j];
    M(4,0) = _x[i];
    M(4,1) = _y[j-1];
    M(5,0) = _x[i+1];
    M(5,1) = _y[j-1];

    M_ind(0,0) = i;
    M_ind(0,1) = j+1;
    M_ind(1,0) = i+1;
    M_ind(1,1) = j+1;
    M_ind(2,0) = i;
    M_ind(2,1) = j;
    M_ind(3,0) = i+1;
    M_ind(3,1) = j;
    M_ind(4,0) = i;
    M_ind(4,1) = j-1;
    M_ind(5,0) = i+1;
    M_ind(5,1) = j-1;

    xj = _grid->get_xinter(i,j);
    yj = _y[j];
    _grid->discret_scheme(i,j,fix,fiy,xj,yj,M, M_ind);

    xk =_grid->get_xk(i,j);
    yk =_grid->get_yk(i,j);
    indk = _grid->get_indk(i,j);

    xi =_grid->get_xi(i,j);
    yi =_grid->get_yi(i,j);
    indi = _grid->get_indi(i,j);

    // calcul des coef du gradient (interpolation lineaire sur les 3 points)
    denom = (xj-xk)*(yj-yi)-(xj-xi)*(yj-yk);
    alphaj = (yk-yi)/denom;
    betaj = (xi-xk)/denom;

    denom = (xi-xk)*(yi-yj)-(xi-xj)*(yi-yk);
    alphai = (yk-yj)/denom;
    betai = (xj-xk)/denom;

    denom = (xk-xj)*(yk-yi)-(xk-xi)*(yk-yj);
    alphak = (yj-yi)/denom;
    betak = (xi-xj)/denom;

    // indk
    mat1[m] = _grid->get_tx2(i,j)-1;
    mat2[m] = indk-1;
    mat3[m] = _grid->coef(_grid->get_xinter(i,j),_y[j])*(alphak*fix + betak*fiy);
    m = m +1;

    // indi
    mat1[m] = _grid->get_tx2(i,j)-1;
    mat2[m] = indi-1;
    mat3[m] = _grid->coef(_grid->get_xinter(i,j),_y[j])*(alphai*fix + betai*fiy);
    m = m +1;

    // int(i,i+1),j
    mat1[m] = _grid->get_tx2(i,j)-1;
    mat2[m] = _grid->get_tx2(i,j)-1;
    mat3[m] = ((0.5+_eta)/0.5)*1 + _grid->coef(_grid->get_xinter(i,j),_y[j])*(alphaj*fix + betaj*fiy);
    m = m +1;

    // int(i,i+1),j
    mat1[m] = _grid->get_tx2(i,j)-1;
    mat2[m] = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i,j))-1;
    mat3[m] = -((0.5+_eta)/0.5)*1.;
    m = m +1;

}

void HeatAssemblerforme::electrodefluxy_forme(int i, int j){

  pi = acos(-1.);
 int l,k;
  fix=0;
  fiy=0;

  //calcul normale sur le point d'interface
  fix = (_grid->get_u(i,j)*(_y[j+1]-_grid->get_yinter(i,j))-_grid->get_u(i,j+1)*(_y[j]-_grid->get_yinter(i,j)))/dy;
  fiy = (_grid->get_v(i,j)*(_y[j+1]-_grid->get_yinter(i,j))-_grid->get_v(i,j+1)*(_y[j]-_grid->get_yinter(i,j)))/dy;
  fiix = fix/sqrt(pow(fix,2)+pow(fiy,2));
  fiiy = fiy/sqrt(pow(fix,2)+pow(fiy,2));
  fix = fiix;
  fiy = fiiy;
  
    Eigen:: MatrixXd M, M_ind;
    M.resize(6,2);
    M_ind.resize(6,2);

    M(0,0) = _x[i-1];
    M(0,1) = _y[j+1];
    M(1,0) = _x[i-1];
    M(1,1) = _y[j];
    M(2,0) = _x[i];
    M(2,1) = _y[j+1];
    M(3,0) = _x[i];
    M(3,1) = _y[j];
    M(4,0) = _x[i+1];
    M(4,1) = _y[j+1];
    M(5,0) = _x[i+1];
    M(5,1) = _y[j];

    M_ind(0,0) = i-1;
    M_ind(0,1) = j+1;
    M_ind(1,0) = i-1;
    M_ind(1,1) = j;
    M_ind(2,0) = i;
    M_ind(2,1) = j+1;
    M_ind(3,0) = i;
    M_ind(3,1) = j;
    M_ind(4,0) = i+1;
    M_ind(4,1) = j+1;
    M_ind(5,0) = i+1;
    M_ind(5,1) = j;

    xj = _x[i];
    yj = _grid->get_yinter(i,j);
    _grid->discret_scheme(i,j,fix,fiy,xj,yj,M, M_ind);

    xk =_grid->get_xk(i,j);
    yk =_grid->get_yk(i,j);
    indk = _grid->get_indk(i,j);

    xi =_grid->get_xi(i,j);
    yi =_grid->get_yi(i,j);
    indi = _grid->get_indi(i,j);

    denom = (xj-xk)*(yj-yi)-(xj-xi)*(yj-yk);
    alphaj = (yk-yi)/denom;
    betaj = (xi-xk)/denom;

    denom = (xi-xk)*(yi-yj)-(xi-xj)*(yi-yk);
    alphai = (yk-yj)/denom;
    betai = (xj-xk)/denom;

    denom = (xk-xj)*(yk-yi)-(xk-xi)*(yk-yj);
    alphak = (yj-yi)/denom;
    betak = (xi-xj)/denom;

    // indk
    mat1[m] = _grid->get_ty2(i,j)-1;
    mat2[m] = indk-1;
    mat3[m] = _grid->coef(_x[i],_grid->get_yinter(i,j))*(alphak*fix + betak*fiy);
    m = m +1;

    // indi
    mat1[m] = _grid->get_ty2(i,j)-1;
    mat2[m] = indi-1;
    mat3[m] = _grid->coef(_x[i],_grid->get_yinter(i,j))*(alphai*fix + betai*fiy);
    m = m +1;

    // int(i,i+1),j
    mat1[m] = _grid->get_ty2(i,j)-1;
    mat2[m] = _grid->get_ty2(i,j)-1;
    mat3[m] = ((0.5+_eta)/0.5)*1 + _grid->coef(_x[i],_grid->get_yinter(i,j))*(alphaj*fix + betaj*fiy);
    m = m +1;

    // Um
    mat1[m] = _grid->get_ty2(i,j)-1;
    mat2[m] = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i,j))-1;
    mat3[m] = -((0.5+_eta)/0.5)*1.;
    m = m +1;
}


void HeatAssemblerforme::BuildLaplacianMatrix_forme(VectorXd alpha)
{

  double b;

  // _grid->BuildT(alpha);
  // _grid->inter(alpha);
  // _grid->DLevelset(alpha);
  // _grid->Zonelectrode(alpha);


  _lap_mat.resize(Nx*Ny+_grid->Get_nbsigne()+Ne,Nx*Ny+_grid->Get_nbsigne()+Ne);

  coeff = 0;
  m=0;

  for (int j=0; j<Ny; ++j){
    for (int i=0; i<Nx; ++i){

      mat1[m]= _grid->ind(i,j)-1;
      mat2[m]= _grid->ind(i,j)-1;
      mat3[m]=0;

      coeff = _grid->coef(_x[i]+dx/2,_y[j]);
      if (_grid->get_xix2_2(i,j) != 0)  coeff = _grid->coef(_grid->get_xinter(i,j)/2 + _x[i]/2, _y[j]);
      mat3[m] =  coeff*( -1./(dx*(1. - _grid->get_xix2_2(i,j)) +_grid->get_xix2_2(i,j)*abs(_grid->get_xinter(i,j)-_x[i])))/dx;

      coeff = _grid->coef(_x[i]-dx/2,_y[j]);
      if(i==0) {
        mat3[m] = mat3[m] + coeff*( -1./(dx*(1.)))/dx;
      }else{
        if (_grid->get_xix2_2(i-1,j) != 0)  coeff = _grid->coef(_grid->get_xinter(i-1,j)/2 + _x[i]/2, _y[j]);
        mat3[m] = mat3[m] + coeff*( -1./(dx*(1. -_grid->get_xix2_2(i-1,j)) +_grid->get_xix2_2(i-1,j)*abs(_x[i]- _grid->get_xinter(i-1,j))))/dx;
      }

      coeff = _grid->coef(_x[i],_y[j] + dy/2.);
      if (_grid->get_xiy2_2(i,j) != 0)  coeff = _grid->coef(_x[i], _grid->get_yinter(i,j)/2 + _y[j]/2);
      mat3[m] = mat3[m] + coeff*( -1./(dy*(1. - _grid->get_xiy2_2(i,j)) +_grid->get_xiy2_2(i,j)*abs(_grid->get_yinter(i,j)-_y[j])))/dx;

      coeff = _grid->coef(_x[i],_y[j]-dy/2);
      if (j==0){
        mat3[m] = mat3[m] + coeff*( -1./(dy*(1.)))/dx;
      }else{
        if (_grid->get_xiy2_2(i,j-1) != 0)  coeff = _grid->coef(_x[i], _grid->get_yinter(i,j-1)/2. + _y[j]/2.);
        mat3[m] = mat3[m] + coeff*( -1./(dy*(1. -_grid->get_xiy2_2(i,j-1)) +_grid->get_xiy2_2(i,j-1)*abs(_y[j]- _grid->get_yinter(i,j-1))))/dx;
      }

      m=m+1;

      if (i-1 >= 0){
        coeff = _grid->coef(_x[i]-dx/2.,_y[j]);
        if (_grid->get_xix2_2(i-1,j) != 0)  coeff = _grid->coef(_grid->get_xinter(i-1,j)/2. + _x[i]/2., _y[j]);
        mat1[m] = _grid->ind(i,j)-1;
        mat2[m] = (1-_grid->get_xix2_2(i-1,j))*_grid->ind(i-1,j) + _grid->get_xix2_2(i-1,j)*_grid->get_tx2(i-1,j)-1;
        mat3[m] =  coeff*((1-_grid->get_xix2_2(i-1,j))/(dx*dx)  + _grid->get_xix2_2(i-1,j)/(abs(_x[i]-_grid->get_xinter(i-1,j))*dx));
        m = m +1;
      }
      // i+1,j
      if (i+1 <= Nx-1){
        coeff = _grid->coef(_x[i]+dx/2.,_y[j]);
        if (_grid->get_xix2_2(i,j) != 0)  coeff = _grid->coef(_grid->get_xinter(i,j)/2. + _x[i]/2, _y[j]);
        mat1[m] = _grid->ind(i,j)-1;
        mat2[m] = ((1-  _grid->get_xix2_2(i,j))*_grid->ind(i+1,j) + _grid->get_xix2_2(i,j)*_grid->get_tx2(i,j))-1;
        mat3[m] = coeff*((1-_grid->get_xix2_2(i,j))/(dx*dx) + _grid->get_xix2_2(i,j)/(abs( _grid->get_xinter(i,j)-_x[i])*dx));
        m = m +1;
      }
      //i,j-1
      if (j-1 >= 0){
        coeff = _grid->coef(_x[i],_y[j]-dy/2.);
        if (_grid->get_xiy2_2(i,j-1) != 0)  coeff = _grid->coef(_x[i], _grid->get_yinter(i,j-1)/2. + _y[j]/2.);
        mat1[m] = _grid->ind(i,j)-1;
        mat2[m] =  (1-_grid->get_xiy2_2(i,j-1))*_grid->ind(i,j-1) + _grid->get_xiy2_2(i,j-1)*_grid->get_ty2(i,j-1)-1;
        mat3[m] = coeff*((1-_grid->get_xiy2_2(i,j-1))/(dy*dy)  + _grid->get_xiy2_2(i,j-1)/(abs(_y[j]-_grid->get_yinter(i,j-1))*dy));
        m = m +1;
      }
      //i,j+1
      if (j+1 <= Nx-1){
        coeff = _grid->coef(_x[i],_y[j] +dy/2.);
        if (_grid->get_xiy2_2(i,j) != 0)  coeff = _grid->coef(_x[i], _grid->get_yinter(i,j)/2. + _y[j]/2.);
        mat1[m] = _grid->ind(i,j)-1;
        mat2[m] =  ((1-  _grid->get_xiy2_2(i,j))*_grid->ind(i,j+1) + _grid->get_xiy2_2(i,j)*_grid->get_ty2(i,j))-1;
        mat3[m] = coeff*((1-_grid->get_xiy2_2(i,j))/(dy*dy) + _grid->get_xiy2_2(i,j)/(abs( _grid->get_yinter(i,j)-_y[j])*dy));
        m = m +1;
      }
    }
  }




  for (int j=1; j< Ny-1; ++j){
    for (int i=1; i<Nx-1; ++i){

      // s'il y a un point d'interface entre i et i+1
      if ((_grid->get_tx2(i,j))!= 0) {
        if (_grid->get_zone(_grid->get_tx2(i,j)) == 0) {
          fluxx_forme(i,j);
        }else{
          electrodefluxx_forme(i,j);
        }
      }

      // s'il y a un point d'interface entre i et i+1
      if ((_grid->get_ty2(i,j))!= 0) {
        if (_grid->get_zone(_grid->get_ty2(i,j)) == 0) {
          fluxy_forme(i,j);
        }else{
          electrodefluxy_forme(i,j);
        }
      }


    }
  }



  r_0 = 1;
  // longueur arc de cercle electrode
  pi = acos(-1.);
  somme = 0;


  for(int i=0; i<_grid->getNbPoints(2); i++){
    Um(i) = 0;
  }


  for(int i=1; i<Nx-1; i++){
    for(int j=1; j<Ny-1; j++){

      if (_grid->get_tx2(i,j) != 0) {
        if (_grid->get_zone(_grid->get_tx2(i,j))  != 0) {
          Um(_grid->get_zone(_grid->get_tx2(i,j))-1) =  Um(_grid->get_zone(_grid->get_tx2(i,j))-1) + dx*dy*_grid->Delta(i,j,alpha);
          mat1(m) = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i,j))-1;
          mat2(m) =  _grid->ind(i,j)-1;
          mat3(m) = -dx*dy*_grid->Delta(i,j,alpha);
          m = m +1;
        }
      }else{
        if (_grid->get_tx2(i-1,j) != 0) {
          if (_grid->get_zone(_grid->get_tx2(i-1,j))  != 0) {
            Um(_grid->get_zone(_grid->get_tx2(i-1,j))-1) =  Um(_grid->get_zone(_grid->get_tx2(i-1,j))-1) + dx*dy*_grid->Delta(i,j,alpha);
            mat1(m) = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i-1,j))-1;
            mat2(m) =  _grid->ind(i,j)-1;
            mat3(m) = -dx*dy*_grid->Delta(i,j,alpha);
            m = m +1;
          }
        }else{
          if (_grid->get_ty2(i,j) != 0) {
            if (_grid->get_zone(_grid->get_ty2(i,j))  != 0) {
              Um(_grid->get_zone(_grid->get_ty2(i,j))-1) =  Um(_grid->get_zone(_grid->get_ty2(i,j))-1) + dx*dy*_grid->Delta(i,j,alpha);
              mat1(m) = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i,j))-1;
              mat2(m) =  _grid->ind(i,j)-1;
              mat3(m) = -dx*dy*_grid->Delta(i,j,alpha);
              m = m +1;
            }
          }else{
            if (_grid->get_ty2(i,j-1) != 0) {
              if (_grid->get_zone(_grid->get_ty2(i,j-1))  != 0) {
                Um(_grid->get_zone(_grid->get_ty2(i,j-1))-1) =  Um(_grid->get_zone(_grid->get_ty2(i,j-1))-1) + dx*dy*_grid->Delta(i,j,alpha);
                mat1(m) = _grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i,j-1))-1;
                mat2(m) =  _grid->ind(i,j)-1;
                mat3(m) = -dx*dy*_grid->Delta(i,j,alpha);
                m = m +1;
              }
            }
          }
        }
      }

    }
  }

  for(int i=0; i<Ne; i++){
    mat1(m) = _grid->Get_nbsigne() + i;
    mat2(m) = _grid->Get_nbsigne() + i;
    mat3[m] = Um(i);
    //cout<<Um(i)<<endl;
    if (i==0) {mat3[m] = mat3[m] + 0.0000000001;}

    m = m +1;
  }

  m = m-1;



  matt1.resize(m+1);
  matt2.resize(m+1);
  matt3.resize(m+1);

  for(int i=0; i< m+1; i++){
    matt1[i]=mat1[i];
    matt2[i]=mat2[i];
    matt3[i]=mat3[i];
    trp.push_back(Trip(matt1[i],matt2[i],matt3[i]));
  }

  _lap_mat.setFromTriplets(trp.begin(), trp.end());

}

void HeatAssemblerforme::CLEAR()
{
trp.clear();
matt1.resize(0);
matt2.resize(0);
matt3.resize(0);
}

void HeatAssemblerforme::BuildSourceTerm_forme(VectorXd& alpha, Eigen::MatrixXd& u, Eigen::VectorXd& Imm)
{
  _source_term.resize(Nx*Ny+_grid->Get_nbsigne()+Ne);

  for(int i=0; i<Nx*Ny+_grid->Get_nbsigne()+Ne; ++i)
  {
    _source_term[i]=0;
  }

  for(int i=0; i<Nx; ++i)
  {
    for(int j=0; j<Ny; ++j)
    {
      _source_term[_grid->ind(i,j)-1] = 0;

      if (_grid->get_tx2(i,j) != 0) {
        if (_grid->get_zone(_grid->get_tx2(i,j))  == 0) {         
          _source_term[_grid->get_tx2(i,j)-1] = 0;
        }
        else{         
          _source_term[_grid->get_tx2(i,j)-1] = (_eta/0.5)*(u(0,_grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i,j))-1) - u(0,_grid->Get_nbsigne() + _grid->get_tx2(i,j)-1));
        }
      }
      if (_grid->get_ty2(i,j) != 0) {
        if (_grid->get_zone(_grid->get_ty2(i,j))  == 0) {         
          _source_term[_grid->get_ty2(i,j)-1] =  0; // ajouter les coef dee sigma 
        }
        else{
          _source_term[_grid->get_ty2(i,j)-1] = (_eta/0.5)*(u(0,_grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i,j))-1) - u(0,_grid->Get_nbsigne() + _grid->get_ty2(i,j)-1));
        }
      }

    }
  }
  //Imm.resize(Ne);

  for(int i=0; i<Ne; ++i)
  {
    _source_term[_grid->Get_nbsigne() +i+1-1]= -(_eta/(0.5+_eta))*Imm(i);
    //_source_term[_grid->Get_nbsigne() +i+1-1]= +(0.5/(0.5+_eta))*Imm(i);
  }
}

void HeatAssemblerforme::BuildSourceTerm_forme_source(Eigen::VectorXd& alpha, Eigen::MatrixXd& u)
{
  _source_term2.resize(Nx*Ny+_grid->Get_nbsigne()+Ne);

  for(int j=0; j<Ny; ++j)
  {
    for(int i=0; i<Nx; ++i)
    {
      _source_term2[_grid->ind(i,j)-1] = _grid->source(_x[i],_y[j]);
    }
  }

  for(int j=0; j<Ny; ++j)
  {
    for(int i=0; i<Nx; ++i)
    {
      if (_grid->get_tx2(i,j) != 0) {
        if (_grid->get_zone(_grid->get_tx2(i,j))  == 0) {
          fix = (_grid->get_u(i,j)*(_x[i+1]-_grid->get_xinter(i,j))-_grid->get_u(i+1,j)*(_x[i]-_grid->get_xinter(i,j)))/dx;
          fiy = (_grid->get_v(i,j)*(_x[i+1]-_grid->get_xinter(i,j))-_grid->get_v(i+1,j)*(_x[i]-_grid->get_xinter(i,j)))/dx;
          fiix = fix/sqrt(pow(fix,2)+pow(fiy,2));
          fiiy = fiy/sqrt(pow(fix,2)+pow(fiy,2));
          fix = fiix;
          fiy = fiiy;
          _source_term2[_grid->get_tx2(i,j)-1] =_grid->val_du(_grid->get_xinter(i,j),_y[j],fix,fiy);
        }
        else{
          fix = (_grid->get_u(i,j)*(_x[i+1]-_grid->get_xinter(i,j))-_grid->get_u(i+1,j)*(_x[i]-_grid->get_xinter(i,j)))/dx;
          fiy = (_grid->get_v(i,j)*(_x[i+1]-_grid->get_xinter(i,j))-_grid->get_v(i+1,j)*(_x[i]-_grid->get_xinter(i,j)))/dx;
          fiix = fix/sqrt(pow(fix,2)+pow(fiy,2));
          fiiy = fiy/sqrt(pow(fix,2)+pow(fiy,2));
          fix = fiix;
          fiy = fiiy;
          _source_term2[_grid->get_tx2(i,j)-1] = _grid->cond_inter(_grid->get_xinter(i,j),_y[j]) + _grid->val_du(_grid->get_xinter(i,j),_y[j],fix,fiy) + (_eta/0.5)*(u(0,_grid->Get_nbsigne() + _grid->get_zone(_grid->get_tx2(i,j))-1) - u(0,_grid->Get_nbsigne() + _grid->get_tx2(i,j)-1));
        }
      }
      if (_grid->get_ty2(i,j) != 0) {
        if (_grid->get_zone(_grid->get_ty2(i,j))  == 0) {
          fix = (_grid->get_u(i,j)*(_y[j+1]-_grid->get_yinter(i,j))-_grid->get_u(i,j+1)*(_y[j]-_grid->get_yinter(i,j)))/dy;
          fiy = (_grid->get_v(i,j)*(_y[j+1]-_grid->get_yinter(i,j))-_grid->get_v(i,j+1)*(_y[j]-_grid->get_yinter(i,j)))/dy;
          fiix = fix/sqrt(pow(fix,2)+pow(fiy,2));
          fiiy = fiy/sqrt(pow(fix,2)+pow(fiy,2));
          fix = fiix;
          fiy = fiiy;
          _source_term2[_grid->get_ty2(i,j)-1] =  _grid->val_du(_x[i],_grid->get_yinter(i,j),fix,fiy); // ajouter les coef dee sigma 
        }
        else{
          fix = (_grid->get_u(i,j)*(_y[j+1]-_grid->get_yinter(i,j))-_grid->get_u(i,j+1)*(_y[j]-_grid->get_yinter(i,j)))/dy;
          fiy = (_grid->get_v(i,j)*(_y[j+1]-_grid->get_yinter(i,j))-_grid->get_v(i,j+1)*(_y[j]-_grid->get_yinter(i,j)))/dy;
          fiix = fix/sqrt(pow(fix,2)+pow(fiy,2));
          fiiy = fiy/sqrt(pow(fix,2)+pow(fiy,2));
          fix = fiix;
          fiy = fiiy;
          _source_term2[_grid->get_ty2(i,j)-1] = _grid->cond_inter(_x[i],_grid->get_yinter(i,j)) + _grid->val_du(_x[i],_grid->get_yinter(i,j),fix,fiy) + (_eta/0.5)*(u(0,_grid->Get_nbsigne() + _grid->get_zone(_grid->get_ty2(i,j))-1) - u(0,_grid->Get_nbsigne() + _grid->get_ty2(i,j)-1));
        }
      }

    }
  }

 
  Ummm.resize(Ne);
   for(int i=0; i<Ne; i++){
     Ummm(i) = 0;
   }


  for(int i=1; i<Nx-1; i++){
    for(int j=1; j<Ny-1; j++){

      if (_grid->get_tx2(i,j) != 0) {
        if (_grid->get_zone(_grid->get_tx2(i,j))  != 0) {
             Ummm(_grid->get_zone(_grid->get_tx2(i,j))-1) =  Ummm(_grid->get_zone(_grid->get_tx2(i,j))-1) + _grid->cond_inter(_x[i],_y[j])*dx*dy*_grid->Delta(i,j,alpha);
        }
      }else{
        if (_grid->get_tx2(i-1,j) != 0) {
          if (_grid->get_zone(_grid->get_tx2(i-1,j))  != 0) {
             Ummm(_grid->get_zone(_grid->get_tx2(i-1,j))-1) =  Ummm(_grid->get_zone(_grid->get_tx2(i-1,j))-1) + _grid->cond_inter(_x[i],_y[j])*dx*dy*_grid->Delta(i,j,alpha);
          }
        }else{
          if (_grid->get_ty2(i,j) != 0) {
            if (_grid->get_zone(_grid->get_ty2(i,j))  != 0) {
             Ummm(_grid->get_zone(_grid->get_ty2(i,j))-1) =  Ummm(_grid->get_zone(_grid->get_ty2(i,j))-1) + _grid->cond_inter(_x[i],_y[j])*dx*dy*_grid->Delta(i,j,alpha);
            }
          }else{
            if (_grid->get_ty2(i,j-1) != 0) {
              if (_grid->get_zone(_grid->get_ty2(i,j-1))  != 0) {
             Ummm(_grid->get_zone(_grid->get_ty2(i,j-1))-1) =  Ummm(_grid->get_zone(_grid->get_ty2(i,j-1))-1) + _grid->cond_inter(_x[i],_y[j])*dx*dy*_grid->Delta(i,j,alpha);
             }
            }
          }
        }
      }

    }
  }


for(int i=0; i<Ne; ++i){
  _source_term2[_grid->Get_nbsigne() +i+1-1]= 0*Um(i)-Ummm(i);
}



}