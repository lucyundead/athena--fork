//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cooling.cpp
//  \brief Problem generator for testing a user-defined cooling source term 
//
//========================================================================================

// Athena++ headers
#include <cmath>
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Kelvin-Helmholz test
//========================================================================================

static Real den0,vflow1,vflow2,amp,pressure;
static Real kb     = 1.38066e-16;              //erg/K 
static Real mn     = 1.27*1.007*1.660539e-24;  //g
static Real pc     = 3.08567758e18;            //cm

static Real lenu    = pc;                      // * 1 cm
static Real velu    = 1.0e5;                   // * 1 cm/s
static Real rhou    = mn   ;                   // * 1 g/cm^3          
static Real energyu = mn*velu*velu;            // * 1 erg/cm^3, 1erg = 1g*cm^2/s^2
static Real timeu   = lenu/velu;               // * 1 s
static Real pressu  = energyu;                 // * 1 erg/cm^3, (g/(cm*s^2))

static Real GAMMAR = 2.0e-26 ;                 //cooling coefficient

long int iseed = -1;

Real Coolingcoeff(int k, int j, int i, const AthenaArray<Real> &prim)

{

 
 Real nn,LAMBDA,Density,Pressure,Tempture; 
 Real coeff;


               nn = prim(IDN,k,j,i);
          Density = prim(IDN,k,j,i)*rhou;
        Pressure  = prim(IPR,k,j,i)*pressu;
        Tempture  = Pressure*mn/(kb*Density);
 
        //Tempture floor to make sure that no strong cooling 

        if (Tempture <= 10.)
        {

            Tempture = 10.;

        }


        LAMBDA    = GAMMAR*((1.0e7)*exp(-118440./(Tempture+1000.))
                  + 1.4e-2*sqrt(Tempture)*exp(-92./Tempture));
        coeff     = nn*(-GAMMAR+nn*LAMBDA);
    


 return coeff; 

}

void Coolingfunc(MeshBlock *pmb, const Real time, const Real dt,
                 const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
                 AthenaArray<Real> &cons)

{


 Real dt_cgs  = dt * timeu;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
  for (int i=pmb->is; i<=pmb->ie; ++i) {

        cons(IEN,k,j,i) -= Coolingcoeff(k,j,i,prim)*dt_cgs/energyu;

  }}}


return;

}



void Inflow_ix1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)

{


  for(int k=ks; k<=ke; ++k){
  for(int j=js; j<=je; ++j){
    for(int i=1; i<=(NGHOST); ++i){

        prim(IDN,k,j,is-i) = den0;
        prim(IVX,k,j,is-i) = vflow1 + amp*(ran2(&iseed) - 0.5);
        prim(IVY,k,j,is-i) = 0.0    + amp*(ran2(&iseed) - 0.5);
        prim(IVZ,k,j,is-i) = 0.0;
        prim(IPR,k,j,is-i) = pressure;

    }
   }}

return;

}


void Inflow_ox1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)

{


  for(int k=ks; k<=ke; ++k){
  for(int j=js; j<=je; ++j){
    for(int i=1; i<=(NGHOST); ++i){

        prim(IDN,k,j,ie+i) = den0;
        prim(IVX,k,j,ie+i) = -vflow1 + amp*(ran2(&iseed) - 0.5);
        prim(IVY,k,j,ie+i) = 0.0     + amp*(ran2(&iseed) - 0.5);
        prim(IVZ,k,j,ie+i) = 0.0;
        prim(IPR,k,j,ie+i) = pressure;

    }
   }}

return;

}


Real MyTimeStep(MeshBlock *pmb)

{

 Real min_dt=1.0e10;
 Real gm1 = pmb->peos->GetGamma() - 1.0;
 Real energy,press,coolingrate,dt;
 
 for (int k=pmb->ks; k<=pmb->ke; ++k) {
 for (int j=pmb->js; j<=pmb->je; ++j) {
 for (int i=pmb->is; i<=pmb->ie; ++i) {
 
      Real dt;
      press = pmb->phydro->w(IPR,k,j,i);
      energy = press/gm1;
      coolingrate = Coolingcoeff(k,j,i,pmb->phydro->w)*timeu/energyu;
      dt = 0.1*energy/fabs(coolingrate); 
      min_dt = std::min(min_dt, dt);
      
    
 }}}


  return min_dt;

}


void Mesh::InitUserMeshData(ParameterInput *pin)

{

  vflow1 = pin->GetReal("problem","vflow1");
  vflow2 = pin->GetReal("problem","vflow2");
  amp = pin->GetReal("problem","amp");
  pressure = pin->GetReal("problem","press")*kb/pressu;
  den0 = pin->GetReal("problem","den0");

  EnrollUserExplicitSourceFunction(Coolingfunc);
  EnrollUserTimeStepFunction(MyTimeStep);
  EnrollUserBoundaryFunction(INNER_X1, Inflow_ix1);
  EnrollUserBoundaryFunction(OUTER_X1, Inflow_ox1);
  

  return;

}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  int  nx1 = ie - is + 1;
  Real gm1 = peos->GetGamma() - 1.0;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {

    phydro->u(IDN,k,j,i) = den0;

    if ( i <= nx1/2-1) {

       phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*
                              (vflow1 + amp*(ran2(&iseed) - 0.5));

    }
    
    else {

       phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*
                              (-vflow1 + amp*(ran2(&iseed) - 0.5));

    }
    
    phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*
                           amp*(ran2(&iseed) - 0.5);
    phydro->u(IM3,k,j,i) = 0.0;

    phydro->u(IEN,k,j,i) = 0.5*(pow(phydro->u(IM1,k,j,i),2) + pow(phydro->u(IM2,k,j,i),2))
                         /(phydro->u(IDN,k,j,i)) + pressure/gm1;
                         

  }}}


  return;
}
