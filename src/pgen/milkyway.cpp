//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file milkyway.cpp
//  \brief Problem generator for testing a read-in rotating MW potential 
//                               and iso/adia+cooling eos
//                               only works for cartesian coordinate
//========================================================================================

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstring>    // memset
#include <ctime>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../bvals/bvals.hpp"
#include "../utils/utils.hpp"

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//========================================================================================


//================================================
//---global paramters and units in pgen-----------
//================================================

// read-in file properties

#define filename "barforce.dat"

// physical constant

static Real kb     = 1.38066e-16;              //boltzman constant [ergs/K]
static Real Rg     = 8.314e7;                  //gas constant      [ergs/mole/K]
static Real massu  = 1.660538921e-24;          //atomic mass unit  [g]
static Real massp  = 1.67262158e-24;           //mass of proton    [g] 
static Real kpc    = 3.08567758e21;            //kilo parsecs      [cm]
static Real km     = 1.0e5;                    //kilo meters       [cm]
static Real solarm = 1.9891e33;                //mass of the sun   [g]
static Real years  = 3.1556926e7;              //years in seconds  [s]
static Real Gyr    = 3.1556926e16;             //giga year         [s] 
static Real second = 1.0;                      //second            [s]
static Real gconst = 4.302;                    //graviy constant   [kpc/(1.0e6Msun)*(km/s)^2]

// set units 

static Real unitL   = kpc;                      // length unit in       [cm]
static Real unitV   = km/second;                // velocity unit in     [cm/s]
static Real unitM   = 1.0e6*solarm;             // mass unit in         [g]          

static Real unitD   = unitM/pow(unitL,3);       // volume density unit  [g/cm^3]
static Real unitS   = unitM/pow(unitL,2);       // surface density unit [g/cm^2]
static Real unitT   = unitL/unitV;              // time unit            [s]
static Real unitE   = SQR(unitV);               // energy unit          [cm^2/s^2]
static Real unitP   = unitE*unitD;              // pressure unit        [erg/cm^3] , 1erg = 1g*cm^2/s^2


// numerical floors 

static Real dfloor = 1.0e-6;                  // density floor     [code unit unitD]
static Real pfloor = 1.0e-6;                  // pressure floor    [code unit unitP]
static Real vfloor = 2000.  ;                  // velocity floor    [code unit unitV]
static Real tfloor = 10.    ;                  // temperature floor [K]


// problem generator read-in values

static Real den0,rgas,zgas,pressure,initT,phibar,phisp;
static Real potxmax,potymax,potzmax,potdx,potdy,potdz;
static int  potnx,potny,potnz;
static Real rcut, FullBt;

static Real m_sp,pitch_sp,zeta_sp,ri_sp,gamma_sp1,gamma_sp2;
static Real sigma_sp,epsilon_sp,gua_r_sp,gua_rsig_sp;

// problem generator parameters 

static Real GAMMAR = 2.0e-26 ;                 // cooling coefficient
static Real gmw    = 1.27273 ;                 // mean molecular weight
static Real C_cool = 0.5     ;                 // cooling timestep controller

// potential/force arrays 

AthenaArray<float> StaticForX,    StaticForY,    StaticForZ;               // Initial force
AthenaArray<float> StaticMODForX, StaticMODForY, StaticMODForZ;            // Final   force


// user defined functions

static void grav_for_xyz (Real &fz,  Real &fy,  Real &fx, 
                          Real z, Real y, Real x, Real time);
static void grav_for_xyzh(Real &fzl, Real &fz, Real &fzr, 
                          Real &fyl, Real &fy, Real &fyr, 
                          Real &fxl, Real &fx, Real &fxr,
                          Real deltaz, Real deltay, Real deltax, 
                          Real z, Real y, Real x, Real time);

static void spiral_pot(Real &pot, Real &fz, Real &fy, Real &fx,
                       Real z, Real y, Real x, Real time);
static Real Junqueirapot(Real z, Real y, Real x);

static Real LinInterp(const AthenaArray<float> &array, Real z, Real y, Real x);
static Real LinInterp_sym(const AthenaArray<float> &array, Real z, Real y, Real x);

static int sgn(Real v);

//-------sorcue terms-----------

Real Coolingcoeff(int k, int j, int i, const AthenaArray<Real> &prim)

{

 // the cooling function in Inoue & Inutsuka (2008)
 
 Real nn,LAMBDA,Pressure,Tempture; 
 Real coeff;

               nn = prim(IDN,k,j,i)*unitD/(gmw*massp);
        Pressure  = prim(IPR,k,j,i)*unitP;
        Tempture  = Pressure/(kb*nn);
 

        LAMBDA    = GAMMAR*((1.0e7)*exp(-118440./(Tempture+1000.))
                   +1.4e-2*sqrt(Tempture)*exp(-92./Tempture));
        coeff     = nn*(-GAMMAR+nn*LAMBDA);
    
        coeff     = coeff/(unitE*unitD/unitT);


 return coeff; 

}



Real MyTimeStep(MeshBlock *pmb)

{


 Real min_dt=1.0e10;
 Real dt = 1.0e-6;

 if (NON_BAROTROPIC_EOS)
 {

   Real gm1 = pmb->peos->GetGamma() - 1.0;
   Real energy,press,coolingrate;
 
   for (int k=pmb->ks; k<=pmb->ke; ++k) {
   for (int j=pmb->js; j<=pmb->je; ++j) {
   for (int i=pmb->is; i<=pmb->ie; ++i) {
 
      press = pmb->phydro->w(IPR,k,j,i);
      energy = press/gm1;
      coolingrate = Coolingcoeff(k,j,i,pmb->phydro->w);
      dt = C_cool*energy/fabs(coolingrate); 
      min_dt = std::min(min_dt, dt);
    
   }}}
  }

  min_dt = std::min(min_dt, dt);

  return min_dt;

}

void gravacc(MeshBlock *pmb, const Real time, 
     const Real dt,const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bcc, AthenaArray<Real> &cons) {

     Real gm1 = pmb->peos->GetGamma() - 1.0;

     Real x,y,z,rad,deltax,deltay,deltaz;
     Real pot,fx,fy,fz;
     Real potxl,potxr,potyl,potyr,potzl,potzr;
     Real fxl,fxr,fyl,fyr,fzl,fzr;
     Real srcM1,srcM2,srcM3,srcEG,srcEC;

     for (int k=pmb->ks; k<=pmb->ke; ++k) {
     for (int j=pmb->js; j<=pmb->je; ++j) {
     for (int i=pmb->is; i<=pmb->ie; ++i) {

         x  = pmb->pcoord->x1v(i);
         y  = pmb->pcoord->x2v(j);
         z  = pmb->pcoord->x3v(k);
         rad = sqrt(x*x + y*y);

         if (NON_BAROTROPIC_EOS)
         {

             deltax = pmb->pcoord->dx1f(i);
             deltay = pmb->pcoord->dx2f(j);
             deltaz = pmb->pcoord->dx3f(k);

             grav_for_xyzh(fzl,fz,fzr,fyl,fy,fyr,fxl,fx,fxr,
                           deltaz,deltay,deltax,z,y,x,time);

         } else
         {
             grav_for_xyz(fz,fy,fx,z,y,x,time);
         }

         // Momentum terms

         srcM1 = dt*prim(IDN,k,j,i)*fx;
         srcM2 = dt*prim(IDN,k,j,i)*fy;

         cons(IM1,k,j,i) += srcM1;
         cons(IM2,k,j,i) += srcM2;

         if (pmb->block_size.nx3 > 1)
         {
            srcM3 = dt*prim(IDN,k,j,i)*fz;
            cons(IM3,k,j,i) += srcM3;
         }

         // Energy terms 

         if (NON_BAROTROPIC_EOS) 
         {


             srcEG = pmb->phydro->flux[X1DIR](IDN,k,j,i  )*fxl
                    +pmb->phydro->flux[X1DIR](IDN,k,j,i+1)*fxr
                    +pmb->phydro->flux[X2DIR](IDN,k,j  ,i)*fyl
                    +pmb->phydro->flux[X2DIR](IDN,k,j+1,i)*fyr;

             if (pmb->block_size.nx3 > 1)

             {

                srcEG += pmb->phydro->flux[X3DIR](IDN,k  ,j,i)*fzl
                        +pmb->phydro->flux[X3DIR](IDN,k+1,j,i)*fzr;
             }

             srcEC = -Coolingcoeff(k,j,i,prim)*dt;

             cons(IEN,k,j,i) += (0.5*dt*srcEG+srcEC);


         }
  

         // Clean the bounday and set parameter floors


         cons(IDN,k,j,i) = (cons(IDN,k,j,i) > dfloor) ?  cons(IDN,k,j,i) : dfloor;

         cons(IM1,k,j,i) = (cons(IM1,k,j,i) <  vfloor*cons(IDN,k,j,i)) ? cons(IM1,k,j,i) :  vfloor*cons(IDN,k,j,i);
         cons(IM1,k,j,i) = (cons(IM1,k,j,i) > -vfloor*cons(IDN,k,j,i)) ? cons(IM1,k,j,i) : -vfloor*cons(IDN,k,j,i);
         cons(IM2,k,j,i) = (cons(IM2,k,j,i) <  vfloor*cons(IDN,k,j,i)) ? cons(IM2,k,j,i) :  vfloor*cons(IDN,k,j,i);
         cons(IM2,k,j,i) = (cons(IM2,k,j,i) > -vfloor*cons(IDN,k,j,i)) ? cons(IM2,k,j,i) : -vfloor*cons(IDN,k,j,i);
         cons(IM3,k,j,i) = (cons(IM3,k,j,i) <  vfloor*cons(IDN,k,j,i)) ? cons(IM3,k,j,i) :  vfloor*cons(IDN,k,j,i);
         cons(IM3,k,j,i) = (cons(IM3,k,j,i) > -vfloor*cons(IDN,k,j,i)) ? cons(IM3,k,j,i) : -vfloor*cons(IDN,k,j,i);

         if (NON_BAROTROPIC_EOS)
         {

            Real Press_tmp = gm1*(cons(IEN,k,j,i) - 0.5*(SQR(cons(IM1,k,j,i))
                                  +SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)))
                                  /cons(IDN,k,j,i));

            Real Temp  = Press_tmp*unitP/(kb*cons(IDN,k,j,i)*unitD/(gmw*massp));
            Temp = (Temp > tfloor) ?  Temp : tfloor;
            
            Real w_p = (Temp*kb*cons(IDN,k,j,i)*unitD/(gmw*massp))/unitP;
            cons(IEN,k,j,i) = w_p/gm1 + 0.5*(SQR(cons(IM1,k,j,i)) 
                             +SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)))/cons(IDN,k,j,i);

         }


         if (rad > rcut)
         {


            cons(IDN,k,j,i) = pmb->ruser_meshblock_data[0](k,j,i);
            cons(IM1,k,j,i) = pmb->ruser_meshblock_data[1](k,j,i)
                             *pmb->ruser_meshblock_data[0](k,j,i);
            cons(IM2,k,j,i) = pmb->ruser_meshblock_data[2](k,j,i)
                             *pmb->ruser_meshblock_data[0](k,j,i);
            cons(IM3,k,j,i) = pmb->ruser_meshblock_data[3](k,j,i)
                             *pmb->ruser_meshblock_data[0](k,j,i);

            if (NON_BAROTROPIC_EOS) 
            {

               cons(IEN,k,j,i) = pmb->ruser_meshblock_data[4](k,j,i)/gm1
                                +0.5*(SQR(cons(IM1,k,j,i))
                                +SQR(cons(IM2,k,j,i)) + SQR(cons(IM3,k,j,i)))/cons(IDN,k,j,i);
            }

         }

     }}}

     return;

}


void Isolate_ix3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,int ngh)

{

  Real x,y,z,rad;

  for(int k=1 ; k<=ngh; ++k){
  for(int j=js; j<=je ; ++j){
  for(int i=is; i<=ie ; ++i){

//        x = pco->x1v(i);
//        y = pco->x2v(j);
//        z = pco->x3v(k);
//        rad = sqrt(x*x + y*y);

        prim(IDN,ks-k,j,i) =       prim(IDN,ks,j,i);
        prim(IVX,ks-k,j,i) =       prim(IVX,ks,j,i);
        prim(IVY,ks-k,j,i) =       prim(IVY,ks,j,i);

        if (prim(IVZ,ks  ,j,i) >= 0.0) {
            prim(IVZ,ks-k,j,i)  = -prim(IVZ,ks,j,i);
        } else {
            prim(IVZ,ks-k,j,i)  =  prim(IVZ,ks,j,i);
        }

        if (NON_BAROTROPIC_EOS)
            prim(IPR,ks-k,j,i) =  prim(IPR,ks,j,i);

    
  }}}

return;

}


void Isolate_ox3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,int ngh)

{

  Real x,y,z,rad;

  for(int k=1 ; k<=ngh; ++k){
  for(int j=js; j<=je ; ++j){
  for(int i=is; i<=ie ; ++i){

//        x = pco->x1v(i);
//        y = pco->x2v(j);
//        z = pco->x3v(k);
//        rad = sqrt(x*x + y*y);

        prim(IDN,ke+k,j,i) =      prim(IDN,ke,j,i);
        prim(IVX,ke+k,j,i) =      prim(IVX,ke,j,i);
        prim(IVY,ke+k,j,i) =      prim(IVY,ke,j,i);
       
        if (prim(IVZ,ke  ,j,i) <= 0.0) {
            prim(IVZ,ke+k,j,i)  = -prim(IVZ,ke,j,i);
        } else {
            prim(IVZ,ke+k,j,i)  =  prim(IVZ,ke,j,i);
        }

        if (NON_BAROTROPIC_EOS)
           prim(IPR,ke+k,j,i) =   prim(IPR,ke,j,i);

    
  }}}

return;

}


//---------Problem Generator---------

void Mesh::InitUserMeshData(ParameterInput *pin)

{

  //-----read in parameters from file

  den0 = pin->GetReal("problem","den0");
  rgas = pin->GetReal("problem","rgas");
  zgas = pin->GetReal("problem","zgas");
  initT = pin->GetReal("problem","initT");
  phibar = pin->GetReal("problem","phibar");
  phisp  = pin->GetReal("problem","phisp");
  rcut = pin->GetReal("problem","rcut");
  FullBt = pin->GetReal("problem","FullBt");

  potxmax = pin->GetReal("problem","potxmax");
  potymax = pin->GetReal("problem","potymax");
  potzmax = pin->GetReal("problem","potzmax");
  potnx   = pin->GetInteger("problem","potnx");
  potny   = pin->GetInteger("problem","potny");
  potnz   = pin->GetInteger("problem","potnz");

  potdx   = 2.*potxmax/(potnx-1);
  potdy   = 2.*potymax/(potny-1);
  if (potnz > 1) potdz = 2.*potzmax/(potnz-1);


  //-----read-in potential allocation

//  if (Globals::my_rank==0)  // only the master process reads the file 
//  {

    FILE* pFile;
  
    StaticForX.NewAthenaArray(potnz,potny,potnx);
    StaticForY.NewAthenaArray(potnz,potny,potnx);
    StaticForZ.NewAthenaArray(potnz,potny,potnx);
    StaticMODForX.NewAthenaArray(potnz,potny,potnx);
    StaticMODForY.NewAthenaArray(potnz,potny,potnx);
    StaticMODForZ.NewAthenaArray(potnz,potny,potnx);

    float FXI[potnz][potny][potnx];
    float FYI[potnz][potny][potnx];
    float FXF[potnz][potny][potnx];
    float FYF[potnz][potny][potnx];
    float FZI[potnz][potny][potnx];
    float FZF[potnz][potny][potnx];

    if ((pFile = fopen(filename, "rb")) == NULL)
    {
       fprintf(stdout, "Error: Can't open file !\n");
       exit(-2);
    }

    pFile = fopen(filename, "rb");

    fread(FXI, sizeof(FXI), 1, pFile);
    fread(FYI, sizeof(FYI), 1, pFile);
    if (potnz > 1) fread(FZI, sizeof(FZI), 1, pFile);

    fread(FXF, sizeof(FXF), 1, pFile);
    fread(FYF, sizeof(FYF), 1, pFile);
    if (potnz > 1) fread(FZF, sizeof(FZF), 1, pFile);

    fclose(pFile);
  

    for (int k=0; k<potnz; k++){
    for (int j=0; j<potny; j++){
    for (int i=0; i<potnx; i++){

        StaticForX(k,j,i) = FXI[k][j][i];
        StaticForY(k,j,i) = FYI[k][j][i];
        StaticForZ(k,j,i) = FZI[k][j][i];
        StaticMODForX(k,j,i) = FXF[k][j][i];
        StaticMODForY(k,j,i) = FYF[k][j][i];
        StaticMODForZ(k,j,i) = FZF[k][j][i];

    }}}


//  } 

  //----- initial calculations

  phibar = phibar * 1.02269032; //km/s/kpc to Gyr^-1
  phisp  = phisp * 1.02269032; //km/s/kpc to Gyr^-1

//  FullBt = FullBt * 2*PI/phibar/(unitT/Gyr);
  FullBt = 0.1;

  m_sp = 2;
  pitch_sp = 12.5;     // 12.5
  zeta_sp = 800.0;     // 800.
  ri_sp = 8.0;
  gamma_sp1 = 139.5;
  gamma_sp2 = 69.75;
  sigma_sp = 2.35;
  epsilon_sp = 3.8;
  gua_r_sp = 9.0;     // 9.0
  gua_rsig_sp = 1.5;  // 1.5

  //-----enroll physical source terms

  EnrollUserExplicitSourceFunction(gravacc);


  //-----enroll User time step function

//  if (NON_BAROTROPIC_EOS)
//  {  
//
//      EnrollUserTimeStepFunction(MyTimeStep);
//
//  }


  //------enroll user bondary conditions
  
  EnrollUserBoundaryFunction(INNER_X3, Isolate_ix3);
  EnrollUserBoundaryFunction(OUTER_X3, Isolate_ox3);

  return;

}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)

{

    //-----initialize user-defined data field

      int NX1 = (ie-is+1)+2*NGHOST;
      int NX2 = (je-js+1)+2*NGHOST;
      int NX3 = (ke-ks+1)+2*NGHOST;

      AllocateRealUserMeshBlockDataField(5);

      ruser_meshblock_data[0].NewAthenaArray(NX3,NX2,NX1);  //  Initial density
      ruser_meshblock_data[1].NewAthenaArray(NX3,NX2,NX1);  //  Initial Vx
      ruser_meshblock_data[2].NewAthenaArray(NX3,NX2,NX1);  //  Initial Vy
      ruser_meshblock_data[3].NewAthenaArray(NX3,NX2,NX1);  //  Initial Vz
      ruser_meshblock_data[4].NewAthenaArray(NX3,NX2,NX1);  //  Initial Pressure

      return;

}



void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  Real gm1 = peos->GetGamma() - 1.0;
  Real x,y,z,rad,radius;
  Real pot,fx,fy,fz,angvel;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {

    x = pcoord->x1v(i);
    y = pcoord->x2v(j);
    z = pcoord->x3v(k); 

    rad = sqrt(x*x + y*y);
 
    grav_for_xyz(fz,fy,fx,
                 z,y,x,0.0);

    angvel= sqrt( sqrt(fx*fx+fy*fy) *rad );

    if (potnz > 1)   // 3D case
    {
      zgas = 0.018*rad + 0.03;
      phydro->u(IDN,k,j,i) = den0/(2*zgas)*exp(-rad/rgas)*1/SQR(cosh(z/zgas));
    } else
    {
      phydro->u(IDN,k,j,i) = den0*exp(-rad/rgas);
    }

    phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*(-angvel*y/rad);
    
    phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*( angvel*x/rad);

    phydro->u(IM3,k,j,i) = 0.0;

    if (NON_BAROTROPIC_EOS) {

        pressure = phydro->u(IDN,k,j,i)*unitD*kb*initT/(gmw*massp)/unitP;
 

        phydro->u(IEN,k,j,i) = 0.5*(pow(phydro->u(IM1,k,j,i),2) 
                                  + pow(phydro->u(IM2,k,j,i),2)
                                  + pow(phydro->u(IM3,k,j,i),2))
                              /(phydro->u(IDN,k,j,i)) + pressure/gm1;

    }

    // store the initial values
   
 
    ruser_meshblock_data[0](k,j,i) = phydro->u(IDN,k,j,i);
    ruser_meshblock_data[1](k,j,i) = phydro->u(IM1,k,j,i)/phydro->u(IDN,k,j,i);
    ruser_meshblock_data[2](k,j,i) = phydro->u(IM2,k,j,i)/phydro->u(IDN,k,j,i);
    ruser_meshblock_data[3](k,j,i) = phydro->u(IM3,k,j,i)/phydro->u(IDN,k,j,i);


    if (NON_BAROTROPIC_EOS) {
                  
       ruser_meshblock_data[4](k,j,i) = pressure;  

    }


  }}}


  return;
}


void MeshBlock::UserWorkInLoop(void) {

  return;
}


void Mesh::UserWorkAfterLoop(ParameterInput *pin) {

  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {

  return;
}


//------potential functions-------------

static void grav_for_xyz(Real &fz, Real &fy, Real &fx, 
                         Real z, Real y, Real x, Real time)
{

  Real potsp,fzsp,fysp,fxsp;
  Real angle,xrot,yrot,zrot,ti;
  Real fxirot,fxfrot,fyirot,fyfrot;
  Real fxi,fyi,fzi;
  Real fxf,fyf,fzf;
  Real rad2d,rad3d;

  ti = time*(unitT/Gyr);
  angle = -1*ti*phibar;
  rad2d = sqrt(SQR(x)+SQR(y));

  spiral_pot(potsp,fzsp,fysp,fxsp,z,y,x,time);

  if (rad2d <= rcut)
  {
     
     xrot = x*cos(angle)-y*sin(angle);
     yrot = x*sin(angle)+y*cos(angle);
     zrot = z;
  
     fxirot = LinInterp(StaticForX,zrot,yrot,xrot);
     fyirot = LinInterp(StaticForY,zrot,yrot,xrot);
     fzi    = LinInterp(StaticForZ,zrot,yrot,xrot)*sgn(z);

     fxfrot = LinInterp(StaticMODForX,zrot,yrot,xrot);
     fyfrot = LinInterp(StaticMODForY,zrot,yrot,xrot);
     fzf    = LinInterp(StaticMODForZ,zrot,yrot,xrot)*sgn(z);

     fxi = fxirot*cos(-angle)-fyirot*sin(-angle);
     fyi = fxirot*sin(-angle)+fyirot*cos(-angle);
     
     fxf = fxfrot*cos(-angle)-fyfrot*sin(-angle);
     fyf = fxfrot*sin(-angle)+fyfrot*cos(-angle);

  } else
  {

     fxi = LinInterp(StaticForX,z,y,x);
     fyi = LinInterp(StaticForY,z,y,x);    
     fzi = LinInterp(StaticForZ,z,y,x)*sgn(z);

     fxf = LinInterp(StaticMODForX,z,y,x);
     fyf = LinInterp(StaticMODForY,z,y,x);
     fzf = LinInterp(StaticMODForZ,z,y,x)*sgn(z);

  }

  if (time <= FullBt)
  {

      fx = (fxf+fxsp)*time/FullBt
          +(fxi     )*(1.-time/FullBt);
      fy = (fyf+fysp)*time/FullBt
          +(fyi     )*(1.-time/FullBt);
      fz = (fzf+fzsp)*time/FullBt
          +(fzi     )*(1.-time/FullBt);
  } else
  {
      fx = fxf+fxsp;
      fy = fyf+fysp;
      fz = fzf+fzsp;    
  }

  return;

}

static void grav_for_xyzh(Real &fzl, Real &fz, Real &fzr, Real &fyl, Real &fy, Real &fyr,
                          Real &fxl, Real &fx, Real &fxr, Real deltaz, Real deltay, Real deltax,
                          Real z, Real y, Real x, Real time)
{

  Real potsp,fzsp,fysp,fxsp;
  Real potspxl,potspxr,potspyl,potspyr,potspzl,potspzr;
  Real fzspl,fzspr,fyspl,fyspr,fxspl,fxspr;
  Real angle,xrot,yrot,zrot,ti;
  Real xl,xr,yl,yr,zl,zr;
  Real fxi,fyi,fzi,fxirot,fyirot;
  Real fxf,fyf,fzf,fxfrot,fyfrot;
  Real fxil,fxir,fyil,fyir,fzil,fzir;
  Real fxfl,fxfr,fyfl,fyfr,fzfl,fzfr;
  Real rad2d,rad3d,tmp1,tmp2,tmp3;

  ti = time*(unitT/Gyr);
  angle = -1*ti*phibar;
  rad2d = sqrt(SQR(x)+SQR(y));

  xl = x - 0.5*deltax;
  xr = x + 0.5*deltax;
  yl = y - 0.5*deltay;
  yr = y + 0.5*deltay;
  zl = z - 0.5*deltaz;
  zr = z + 0.5*deltaz;

  spiral_pot(potsp  ,fzsp,fysp,fxsp ,z ,y ,x ,time);

  spiral_pot(potspxl,tmp1,tmp2,fxspl,z ,y ,xl,time);
  spiral_pot(potspxr,tmp1,tmp2,fxspr,z ,y ,xr,time);

  spiral_pot(potspyl,tmp1,fyspl,tmp2,z ,yl,x ,time);
  spiral_pot(potspyr,tmp1,fyspr,tmp2,z ,yr,x ,time);

  spiral_pot(potspzl,fzspl,tmp1,tmp2,zl,y ,x ,time);
  spiral_pot(potspzr,fzspr,tmp1,tmp2,zr,y ,x ,time);

  if (rad2d <= rcut)
  {

     xrot = x*cos(angle)-y*sin(angle);
     yrot = x*sin(angle)+y*cos(angle);
     zrot = z;

     fxirot = LinInterp(StaticForX,zrot,yrot,xrot);
     fyirot = LinInterp(StaticForY,zrot,yrot,xrot);
     fzi    = LinInterp(StaticForZ,zrot,yrot,xrot)*sgn(z);

     fxfrot = LinInterp(StaticMODForX,zrot,yrot,xrot);
     fyfrot = LinInterp(StaticMODForY,zrot,yrot,xrot);
     fzf  = LinInterp(StaticMODForZ,zrot,yrot,xrot)*sgn(z);

     fxi = fxirot*cos(-angle)-fyirot*sin(-angle);
     fyi = fxirot*sin(-angle)+fyirot*cos(-angle);

     fxf = fxfrot*cos(-angle)-fyfrot*sin(-angle);
     fyf = fxfrot*sin(-angle)+fyfrot*cos(-angle);

     fzil = LinInterp(StaticForZ,zl,yrot,xrot)*sgn(z);
     fzir = LinInterp(StaticForZ,zr,yrot,xrot)*sgn(z);
     fzfl = LinInterp(StaticMODForZ,zl,yrot,xrot)*sgn(z);
     fzfr = LinInterp(StaticMODForZ,zr,yrot,xrot)*sgn(z);

     xrot = xl*cos(angle)-y*sin(angle);
     yrot = xl*sin(angle)+y*cos(angle);
     fxirot = LinInterp(StaticForX,zrot,yrot,xrot);
     fyirot = LinInterp(StaticForY,zrot,yrot,xrot);
     fxfrot = LinInterp(StaticMODForX,zrot,yrot,xrot);
     fyfrot = LinInterp(StaticMODForY,zrot,yrot,xrot);
     fxil = fxirot*cos(-angle)-fyirot*sin(-angle);
     fxfl = fxfrot*cos(-angle)-fyfrot*sin(-angle);
     
     xrot = xr*cos(angle)-y*sin(angle);
     yrot = xr*sin(angle)+y*cos(angle);
     fxirot = LinInterp(StaticForX,zrot,yrot,xrot);
     fyirot = LinInterp(StaticForY,zrot,yrot,xrot);
     fxfrot = LinInterp(StaticMODForX,zrot,yrot,xrot);
     fyfrot = LinInterp(StaticMODForY,zrot,yrot,xrot);
     fxir = fxirot*cos(-angle)-fyirot*sin(-angle);
     fxfr = fxfrot*cos(-angle)-fyfrot*sin(-angle);

     xrot = x*cos(angle)-yl*sin(angle);
     yrot = x*sin(angle)+yl*cos(angle);
     fxirot = LinInterp(StaticForX,zrot,yrot,xrot);
     fyirot = LinInterp(StaticForY,zrot,yrot,xrot);
     fxfrot = LinInterp(StaticMODForX,zrot,yrot,xrot);
     fyfrot = LinInterp(StaticMODForY,zrot,yrot,xrot);
     fyil = fxirot*sin(-angle)+fyirot*cos(-angle);
     fyfl = fxfrot*sin(-angle)+fyfrot*cos(-angle);

     xrot = x*cos(angle)-yr*sin(angle);
     yrot = x*sin(angle)+yr*cos(angle);
     fxirot = LinInterp(StaticForX,zrot,yrot,xrot);
     fyirot = LinInterp(StaticForY,zrot,yrot,xrot);
     fxfrot = LinInterp(StaticMODForX,zrot,yrot,xrot);
     fyfrot = LinInterp(StaticMODForY,zrot,yrot,xrot);
     fyir = fxirot*sin(-angle)+fyirot*cos(-angle);
     fyfr = fxfrot*sin(-angle)+fyfrot*cos(-angle);

  } else
  {

     fxi = LinInterp(StaticForX,z ,y ,x );
     fyi = LinInterp(StaticForY,z ,y ,x );
     fzi = LinInterp(StaticForZ,z ,y ,x )*sgn(z);
     fxil= LinInterp(StaticForX,z ,y ,xl);
     fxir= LinInterp(StaticForX,z ,y ,xr);
     fyil= LinInterp(StaticForY,z ,yl,x );
     fyir= LinInterp(StaticForY,z ,yr,x );
     fzil= LinInterp(StaticForZ,zl,y ,x )*sgn(z);
     fzir= LinInterp(StaticForZ,zr,y ,x )*sgn(z);

     fxf = LinInterp(StaticMODForX,z ,y ,x );
     fyf = LinInterp(StaticMODForY,z ,y ,x );
     fzf = LinInterp(StaticMODForZ,z ,y ,x )*sgn(z);
     fxfl= LinInterp(StaticMODForX,z ,y ,xl);
     fxfr= LinInterp(StaticMODForX,z ,y ,xr);
     fyfl= LinInterp(StaticMODForY,z ,yl,x );
     fyfr= LinInterp(StaticMODForY,z ,yr,x );
     fzfl= LinInterp(StaticMODForZ,zl,y ,x )*sgn(z);
     fzfr= LinInterp(StaticMODForZ,zr,y ,x )*sgn(z);

  }

  if (time <= FullBt)
  {
      fxl= (fxfl+fxspl)*time/FullBt
          +(fxil      )*(1.-time/FullBt);
      fx = (fxf +fxsp )*time/FullBt
          +(fxi       )*(1.-time/FullBt);
      fxr= (fxfr+fxspr)*time/FullBt
          +(fxir      )*(1.-time/FullBt);
      fyl= (fyfl+fyspl)*time/FullBt
          +(fyil      )*(1.-time/FullBt);
      fy = (fyf +fysp )*time/FullBt
          +(fyi       )*(1.-time/FullBt);
      fyr= (fyfr+fyspr)*time/FullBt
          +(fyir      )*(1.-time/FullBt);
      fzl= (fzfl+fzspl)*time/FullBt
          +(fzil      )*(1.-time/FullBt);
      fz = (fzf +fzsp )*time/FullBt
          +(fzi       )*(1.-time/FullBt);
      fzr= (fzfr+fzspr)*time/FullBt
          +(fzir      )*(1.-time/FullBt);

  } else
  {
      fxl= fxfl+fxspl;
      fx = fxf +fxsp;
      fxr= fxfr+fxspr;
      fyl= fyfl+fyspl;
      fy = fyf +fysp;
      fyr= fyfr+fyspr;
      fzl= fzfl+fzspl;
      fz = fzf +fzsp;
      fzr= fzfr+fzspr;
  }

  return;

}

static void spiral_pot(Real &pot, Real &fz, Real &fy, Real &fx,
                       Real z, Real y, Real x, Real time)
{

  Real angle,xrot,yrot,zrot,delta,ti;
  Real fxsprot,fysprot,fzsprot;

  ti = time*(unitT/Gyr);
  angle = -1*ti*phisp;

  xrot = x*cos(angle)-y*sin(angle);
  yrot = x*sin(angle)+y*cos(angle);
  zrot = z;

//--numerically calculate the force

  delta = 0.001;

  fxsprot = -(Junqueirapot(zrot,yrot,xrot+delta) - 
              Junqueirapot(zrot,yrot,xrot-delta))/(2.*delta);
  fysprot = -(Junqueirapot(zrot,yrot+delta,xrot) - 
              Junqueirapot(zrot,yrot-delta,xrot))/(2.*delta);
  fz      = -(Junqueirapot(zrot+delta,yrot,xrot) - 
              Junqueirapot(zrot-delta,yrot,xrot))/(2.*delta);

  
  fx = fxsprot*cos(-angle)-fysprot*sin(-angle);
  fy = fxsprot*sin(-angle)+fysprot*cos(-angle);

  pot = Junqueirapot(xrot,yrot,zrot);


}


static Real Junqueirapot(Real z, Real y, Real x) 

{

  Real rad,phi_sp,fm1,fm2,k,Junpot;

  rad = sqrt(SQR(x)+SQR(y));

  if ( y >= 0.0) {

     phi_sp = acos(-x/rad);

  }

  else {

     phi_sp = 2.0*PI-acos(-x/rad);

  }


  fm1 = 2.*(1/tan(pitch_sp*PI/180.)*log(rad/ri_sp)+(gamma_sp1*PI/180.));
  fm2 = 2.*(1/tan(pitch_sp*PI/180.)*log(rad/ri_sp)+(gamma_sp2*PI/180.));

  k = 2./(rad*tan(pitch_sp*PI/180.));

  if ( rad >= gua_r_sp ) {

     Junpot = (-1*zeta_sp*rad*exp(-(SQR(rad)/SQR(sigma_sp))
             *(1-cos(m_sp*phi_sp-fm1))-rad/epsilon_sp-fabs(k*z))
               -1*zeta_sp*rad*exp(-(SQR(rad)/SQR(sigma_sp))
             *(1-cos(m_sp*phi_sp-fm2))-rad/epsilon_sp-fabs(k*z)));

  }

  else {


     Junpot = (-1*zeta_sp*rad*exp(-(SQR(rad)/SQR(sigma_sp))
              *(1-cos(m_sp*phi_sp-fm1))-rad/epsilon_sp-fabs(k*z))
              *exp(-SQR(rad-gua_r_sp)/(2.0*SQR(gua_rsig_sp)))
               -1*zeta_sp*rad*exp(-(SQR(rad)/SQR(sigma_sp))
              *(1-cos(m_sp*phi_sp-fm2))-rad/epsilon_sp-fabs(k*z))
              *exp(-SQR(rad-gua_r_sp)/(2.0*SQR(gua_rsig_sp))));

  }


  return Junpot;
}

//------interpolation functions-------------

static Real LinInterp_sym(const AthenaArray<float> &array, Real z, Real y, Real x)
{

  int xindex,yindex,zindex;
  Real c00,c01,c10,c11,c0,c1,xd,yd,zd,result;


  if (potnz > 1)               // 3D case

  {  

         xindex = int(floor((x+potxmax)/potdx));
         yindex = int(floor((y+potymax)/potdy));
         zindex = int(floor((z+potzmax)/potdz));

         xd  = (x - (potdx*(xindex) - potxmax))/potdx;
         yd  = (y - (potdy*(yindex) - potymax))/potdy;
         zd  = (z - (potdz*(zindex) - potzmax))/potdz;

         c00 = array(zindex  ,yindex  ,xindex  )*(1.-xd) 
              +array(zindex  ,yindex  ,xindex+1)*xd;
         c10 = array(zindex  ,yindex+1,xindex  )*(1.-xd) 
              +array(zindex  ,yindex+1,xindex+1)*xd;
         c01 = array(zindex+1,yindex  ,xindex  )*(1.-xd) 
              +array(zindex+1,yindex  ,xindex+1)*xd;
         c11 = array(zindex+1,yindex+1,xindex  )*(1.-xd) 
              +array(zindex+1,yindex+1,xindex+1)*xd;

         c0 = c00*(1.-yd) + c10*yd;
         c1 = c01*(1.-yd) + c11*yd;

         result = c0*(1.-zd) + c1*zd;

  } else                           // 2D case 

  {

         xindex = int(floor((x+potxmax)/potdx));
         yindex = int(floor((y+potymax)/potdy));

         xd  = (x - (potdx*(xindex) - potxmax))/potdx;
         yd  = (y - (potdy*(yindex) - potymax))/potdy;

         c00 = array(0,yindex  ,xindex  )*(1.-xd)
              +array(0,yindex  ,xindex+1)*xd;
         c10 = array(0,yindex+1,xindex  )*(1.-xd)
              +array(0,yindex+1,xindex+1)*xd;

         result = c00*(1.-yd) + c10*yd;
  }


  return result;                                                    

}


// Inteporlate functions that only works for axisymmetric, non-uniform array
// for now I assume z <-> -z symmetric, and z is non-uniform 
// i.e. dz = 10 pc from 0 to 100 pc, and dz = 100pc from 100 to 500 pc


static Real LinInterp(const AthenaArray<float> &array, Real z, Real y, Real x)
{

  int xindex,yindex,zindex;
  Real c00,c01,c10,c11,c0,c1,xd,yd,zd,result;

  Real potdz1 = 0.01;
  Real potdz2 = 0.1 ;
  Real z_critical = 0.1;

  if (potnz > 1)               // 3D case

  {  

         xindex = int(floor((x+potxmax)/potdx));
         yindex = int(floor((y+potymax)/potdy));
         xd  = (x - (potdx*(xindex) - potxmax))/potdx;
         yd  = (y - (potdy*(yindex) - potymax))/potdy;

         if (fabs(z) <= z_critical)
         {
            zindex = int(floor(fabs(z)/potdz1));
            zd  = (fabs(z) - (potdz1*(zindex)))/potdz1;
         }
         else
         {

            zindex = int(floor((fabs(z)-z_critical)/potdz2)) + 10;
            zd  = (fabs(z) - (potdz2*(zindex-10) + z_critical))/potdz2;

         }

         c00 = array(zindex  ,yindex  ,xindex  )*(1.-xd) 
              +array(zindex  ,yindex  ,xindex+1)*xd;
         c10 = array(zindex  ,yindex+1,xindex  )*(1.-xd) 
              +array(zindex  ,yindex+1,xindex+1)*xd;
         c01 = array(zindex+1,yindex  ,xindex  )*(1.-xd) 
              +array(zindex+1,yindex  ,xindex+1)*xd;
         c11 = array(zindex+1,yindex+1,xindex  )*(1.-xd) 
              +array(zindex+1,yindex+1,xindex+1)*xd;

         c0 = c00*(1.-yd) + c10*yd;
         c1 = c01*(1.-yd) + c11*yd;

         result = c0*(1.-zd) + c1*zd;

//         printf("%e,%e,%d,%e,%e,%e\n",z,zd,zindex,c0,c1,result);

  } else                            

  {

         xindex = int(floor((x+potxmax)/potdx));
         yindex = int(floor((y+potymax)/potdy));

         xd  = (x - (potdx*(xindex) - potxmax))/potdx;
         yd  = (y - (potdy*(yindex) - potymax))/potdy;

         c00 = array(0,yindex  ,xindex  )*(1.-xd)
              +array(0,yindex  ,xindex+1)*xd;
         c10 = array(0,yindex+1,xindex  )*(1.-xd)
              +array(0,yindex+1,xindex+1)*xd;

         result = c00*(1.-yd) + c10*yd;

  }


  return result;        

}


int sgn(Real v) 
{
  return (v > 0) - (v < 0);
}                                            
