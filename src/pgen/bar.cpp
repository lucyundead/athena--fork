//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bar.cpp
//  \brief Problem generator for testing a rotating ferrers bar potential 
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

static Real dfloor = 1.0e-4 ;                  // density floor     [code unit unitD]
static Real pfloor = 1.0e-6 ;                  // pressure floor    [code unit unitP]
static Real vfloor = 2000.  ;                  // velocity floor    [code unit unitV]
static Real tfloor = 10.    ;                  // temperature floor [K]


// problem generator read-in values

static Real den0,pressure,initT,phibar;
static Real rcut, FullBt;


// problem generator parameters 

static Real GAMMAR = 2.0e-26 ;                 // cooling coefficient
static Real gmw    = 1.27273 ;                 // mean molecular weight

static Real rb     = 0.3326153155;
static Real rb0    = 0.379;
static Real rhob   = 2.35523767e4;
static Real coefferrs = -60496.875;
static Real aferrs = 5.0;
static Real bferrs = 2.0;

static Real MNmdisk = 2.2 * 1.0e5;
static Real MNadisk = 14.1;
static Real MNbdisk = 0.0;
static Real MHubbbulge  = 3.30276* 1.0e4;
static Real RHubbbulge  = 0.33;
static Real MHubbbulge2  = 4.90678* 1.0e4;
static Real RHubbbulge2  = 0.38;

static Real rHubbMax    = 10.0;
static Real HubbCi  = 1./( log(rHubbMax/RHubbbulge  + sqrt(1.+SQR(rHubbMax/RHubbbulge )))
                     -rHubbMax/(RHubbbulge *sqrt(1.+SQR(rHubbMax/RHubbbulge ))));
static Real HubbCi2 = 1./( log(rHubbMax/RHubbbulge2 + sqrt(1.+SQR(rHubbMax/RHubbbulge2))) 
                     -rHubbMax/(RHubbbulge2*sqrt(1.+SQR(rHubbMax/RHubbbulge2))));

static Real a2,b2,x2,y2,z2,lambda,dlambda;
static Real Fxrot,Fyrot,Fzrot,xrot,yrot,zrot;
static Real sinth, costh, tanth;
static Real w00, w01, w10, w11, w20, w02, w21, w12, w30, w03;
static Real angle;


// user defined functions

static void ferrersbar(Real &pot, Real &fx, Real &fy, Real &fz,
                       Real x, Real y, Real z, Real time);


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
 Real gm1 = pmb->peos->GetGamma() - 1.0;
 Real energy,press,coolingrate,dt;
 
 for (int k=pmb->ks; k<=pmb->ke; ++k) {
 for (int j=pmb->js; j<=pmb->je; ++j) {
 for (int i=pmb->is; i<=pmb->ie; ++i) {
 
      Real dt;
      press = pmb->phydro->w(IPR,k,j,i);
      energy = press/gm1;
      coolingrate = Coolingcoeff(k,j,i,pmb->phydro->w);
      dt = 0.5*energy/fabs(coolingrate); 
      min_dt = std::min(min_dt, dt);
      
    
 }}}


  return min_dt;

}


void gravacc(MeshBlock *pmb, const Real time, 
     const Real dt,const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bcc, AthenaArray<Real> &cons) {

     Real x,y,z,xl,xr,yl,yr,zl,zr;
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

         ferrersbar(pot,fx,fy,fz,
                    x,y,z,time);

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
 
             xl = pmb->pcoord->x1f(i  );
             xr = pmb->pcoord->x1f(i+1);
             yl = pmb->pcoord->x2f(j  );
             yr = pmb->pcoord->x2f(j+1);

             ferrersbar(potxl,fxl,fy,fz,
                        xl,y,z,time);
             ferrersbar(potxr,fxr,fy,fz,
                        xr,y,z,time);
             ferrersbar(potyl,fx,fyl,fz,
                        x,yl,z,time);
             ferrersbar(potyr,fx,fyr,fz,
                        x,yr,z,time);

             srcEG = pmb->phydro->flux[X1DIR](IDN,k,j,i  )*fxl
                    +pmb->phydro->flux[X1DIR](IDN,k,j,i+1)*fxr
                    +pmb->phydro->flux[X2DIR](IDN,k,j  ,i)*fyl
                    +pmb->phydro->flux[X2DIR](IDN,k,j+1,i)*fyr;

             if (pmb->block_size.nx3 > 1)

             {

                zl = pmb->pcoord->x3f(k  );
                zr = pmb->pcoord->x3f(k+1);

                ferrersbar(potzl,fx,fy,fzl,
                           x,y,zl,time);
                ferrersbar(potzr,fx,fy,fzr,
                           x,y,zr,time);


                srcEG += pmb->phydro->flux[X3DIR](IDN,k,j  ,i)*fzl
                        +pmb->phydro->flux[X3DIR](IDN,k+1,j,i)*fzr;
             }

             srcEC = -Coolingcoeff(k,j,i,prim)*dt;

             cons(IEN,k,j,i) += (0.5*dt*srcEG+srcEC);


         }
  


     }}}

     return;

}


void Isolate_ix3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,int ngh)

{


  for(int k=1 ; k<=ngh; ++k){
  for(int j=js; j<=je ; ++j){
  for(int i=is; i<=ie ; ++i){

        prim(IDN,ks-k,j,i) = dfloor;
        prim(IVX,ks-k,j,i) = prim(IVX,ks,j,i);
        prim(IVY,ks-k,j,i) = prim(IVY,ks,j,i);
        prim(IVZ,ks-k,j,i) = 0.0;
        prim(IPR,ks-k,j,i) = pfloor;

    
  }}}

return;

}


void Isolate_ox3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
               Real time, Real dt, int is, int ie, int js, int je, int ks, int ke,int ngh)

{


  for(int k=1 ; k<=ngh; ++k){
  for(int j=js; j<=je ; ++j){
  for(int i=is; i<=ie ; ++i){

        prim(IDN,ke+k,j,i) = dfloor;
        prim(IVX,ke+k,j,i) = prim(IVX,ke,j,i);
        prim(IVY,ke+k,j,i) = prim(IVY,ke,j,i);
        prim(IVZ,ke+k,j,i) = 0.0;
        prim(IPR,ke+k,j,i) = pfloor;

    
  }}}

return;

}


//---------Problem Generator---------

void Mesh::InitUserMeshData(ParameterInput *pin)

{

  //-----read in parameters from file
  
  den0 = pin->GetReal("problem","den0");
  initT = pin->GetReal("problem","initT");
  phibar = pin->GetReal("problem","phibar");
  rcut = pin->GetReal("problem","rcut");
  FullBt = pin->GetReal("problem","FullBt");

  phibar = phibar * 1.02269032; //km/s/kpc to Gyr^-1
  FullBt = FullBt * 2*PI/phibar/(unitT/Gyr);

  //-----enroll physical source terms

  EnrollUserExplicitSourceFunction(gravacc);


  if (NON_BAROTROPIC_EOS)
  {  

     EnrollUserTimeStepFunction(MyTimeStep);

  }


  //------enroll user bondary conditions
  
//  EnrollUserBoundaryFunction(INNER_X3, Isolate_ix3);
//  EnrollUserBoundaryFunction(OUTER_X3, Isolate_ox3);

  return;

}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)

{

    //-----initialize user-defined data field

      int NX1 = (ie-is+1)+2*NGHOST;
      int NX2 = (je-js+1)+2*NGHOST;
      int NX3 = (ke-ks+1)+2*NGHOST;

      AllocateRealUserMeshBlockDataField(9);

      ruser_meshblock_data[0].NewAthenaArray(NX3,NX2,NX1);  //  Initial density
      ruser_meshblock_data[1].NewAthenaArray(NX3,NX2,NX1);  //  Initial Vx
      ruser_meshblock_data[2].NewAthenaArray(NX3,NX2,NX1);  //  Initial Vy
      ruser_meshblock_data[3].NewAthenaArray(NX3,NX2,NX1);  //  Initial Vz
      ruser_meshblock_data[4].NewAthenaArray(NX3,NX2,NX1);  //  Initial Pressure

      ruser_meshblock_data[5].NewAthenaArray(NX3,NX2,NX1);  //  3D potential
      ruser_meshblock_data[6].NewAthenaArray(NX3,NX2,NX1);  //  Initial Vy
      ruser_meshblock_data[7].NewAthenaArray(NX3,NX2,NX1);  //  Initial Vz
      ruser_meshblock_data[8].NewAthenaArray(NX3,NX2,NX1);  //  Initial Pressure



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
 
    ferrersbar(pot,fx,fy,fz,
                   x,y,z,0.0);

    angvel= sqrt( sqrt(fx*fx+fy*fy) *rad );

    phydro->u(IDN,k,j,i) = den0;

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

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {

        Real gam = peos->GetGamma(); 
 
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real rad = sqrt(x*x + y*y);

        Real& u_d  = phydro->u(IDN,k,j,i);
        Real& u_m1 = phydro->u(IM1,k,j,i);
        Real& u_m2 = phydro->u(IM2,k,j,i);
        Real& u_m3 = phydro->u(IM3,k,j,i);

        u_d = (u_d > dfloor) ?  u_d : dfloor;

        u_m1 = (u_m1 <  vfloor*u_d) ?  u_m1 :  vfloor*u_d;
        u_m1 = (u_m1 > -vfloor*u_d) ?  u_m1 : -vfloor*u_d;
        u_m2 = (u_m2 <  vfloor*u_d) ?  u_m2 :  vfloor*u_d;
        u_m2 = (u_m2 > -vfloor*u_d) ?  u_m2 : -vfloor*u_d;
        u_m3 = (u_m3 <  vfloor*u_d) ?  u_m3 :  vfloor*u_d;
        u_m3 = (u_m3 > -vfloor*u_d) ?  u_m3 : -vfloor*u_d;

        if (NON_BAROTROPIC_EOS) {
          Real& w_p  = phydro->w(IPR,k,j,i);
          Real& u_e  = phydro->u(IEN,k,j,i);

          Real Temp  = w_p*unitP/(kb*u_d*unitD/(gmw*massp));
          Temp = (Temp > tfloor) ?  Temp : tfloor;

          w_p = (Temp*kb*u_d*unitD/(gmw*massp))/unitP;

          Real k_e = 0.5*(SQR(u_m1) + SQR(u_m2) + SQR(u_m3))/u_d;
          u_e = w_p/(gam-1.0)+k_e;
        }

        if (rad > rcut) {
        
           u_d  = ruser_meshblock_data[0](k,j,i); 
           u_m1 = ruser_meshblock_data[1](k,j,i)
                 *ruser_meshblock_data[0](k,j,i);
           u_m2 = ruser_meshblock_data[2](k,j,i)
                 *ruser_meshblock_data[0](k,j,i);
           u_m3 = ruser_meshblock_data[3](k,j,i)
                 *ruser_meshblock_data[0](k,j,i);

           if (NON_BAROTROPIC_EOS) {

              Real& u_e  = phydro->u(IEN,k,j,i);
              Real  k_e  = 0.5*(SQR(u_m1)+SQR(u_m2)+SQR(u_m3))/u_d;
              u_e = ruser_meshblock_data[4](k,j,i)/(gam-1.0)+k_e;

           }

        }

    
  }}}

  return;
}



//------potential functions-------------

static void ferrersbar(Real &pot, Real &fx, Real &fy, Real &fz,
                       Real x, Real y, Real z, Real time) 

{

  Real diskterm,diskterm2,rratio,r2term,Hubbforce,rratio2,r2term2,Hubbforce2;
  Real potb,fxb,fyb,fzb,pots,fxs,fys,fzs;

  Real rad = sqrt(x*x+y*y);
  Real radius = sqrt(x*x+y*y+z*z);
  Real dr = 1./radius;
  Real ti = time*(unitT/Gyr);

  angle = -1*ti*phibar;

  xrot = x*cos(angle)-y*sin(angle);
  yrot = x*sin(angle)+y*cos(angle);
  zrot = z;

  x2 = xrot * xrot;
  y2 = yrot * yrot;
  z2 = zrot * zrot;

  a2  = aferrs*aferrs;
  b2  = bferrs*bferrs;


  if ( (y2/a2+x2/b2+z2/b2) > 1.0 ){

     lambda = (-(a2+b2-x2-y2-z2) + sqrt(SQR(a2+b2-x2-y2-z2)
            - 4.0*(a2*b2-y2*b2-x2*a2-z2*a2)))/2.0;
  } 

  else {

     lambda = 0.0;

  }


  dlambda = (b2+lambda)*sqrt(a2+lambda);
  costh = sqrt((b2+lambda)/(a2+lambda));
  sinth = sqrt(1.0-SQR(costh));
  tanth = sinth/costh;


  w00 = log((1.0+sinth)/costh)*2.0/sqrt(a2-b2);
  w10 = 2.0/sqrt((a2-b2)*SQR(a2-b2)) * (log((1.0+sinth)/costh)-sinth);
  w01 = SQR(tanth)*sinth/sqrt(SQR(a2-b2)*(a2-b2)) - w10/2.0;
  w11 = (w01 - w10)/(a2-b2);
  w20 = 2.0/3.0*(1.0/dlambda/(a2+lambda) - w11);
  w02 = 1.0/4.0*(2.0/dlambda/(b2+lambda) - w11);
  w21 = (w11 - w20)/(a2-b2);
  w12 = (w02 - w11)/(a2-b2);
  w30 = 2.0/5.0*(1.0/dlambda/SQR(a2+lambda) - w21);
  w03 = 1.0/6.0*(2.0/dlambda/SQR(b2+lambda) - w12);

  // n=1 case
  potb = coefferrs*(w00+x2*(x2*w02+2.0*y2*w11-2.0*w01)
                       +y2*(y2*w20+2.0*z2*w11-2.0*w10)
                       +z2*(z2*w02+2.0*x2*w02-2.0*w01));

  Fxrot = -4.* coefferrs * xrot * (x2*w02 + y2*w11 + z2*w02 - w01);
  Fyrot = -4.* coefferrs * yrot * (x2*w11 + y2*w20 - w10);
  Fzrot = -4.* coefferrs * zrot * (z2*w02 + y2*w11 + x2*w02 - w01);

  fxb = Fxrot*cos(-angle)-Fyrot*sin(-angle);
  fyb = Fxrot*sin(-angle)+Fyrot*cos(-angle);
  fzb = Fzrot;

  // add a MN disk 

  diskterm = sqrt(z*z+MNbdisk*MNbdisk);
  diskterm2= sqrt( (MNadisk+diskterm)*(MNadisk+diskterm)+rad*rad );


  pot = -gconst*MNmdisk/diskterm2;

  fx  = -x*gconst*MNmdisk/(diskterm2*diskterm2*diskterm2);
  fy  = -y*gconst*MNmdisk/(diskterm2*diskterm2*diskterm2);
  fz  = -z*gconst*MNmdisk*(MNadisk+diskterm)
            /(diskterm2*diskterm2*diskterm2*diskterm);

  // add a Hubble bulge

  rratio    = radius/RHubbbulge;
  r2term    = sqrt(1.0 + rratio*rratio);
  Hubbforce = MHubbbulge*gconst*HubbCi*(dr*dr*log(rratio+r2term)
             - dr/(RHubbbulge*r2term));

  rratio2    = radius/RHubbbulge2;
  r2term2    = sqrt(1.0 + rratio2*rratio2);
  Hubbforce2 = MHubbbulge2*gconst*HubbCi2*(dr*dr*log(rratio2+r2term2) 
             - dr/(RHubbbulge2*r2term2));
  
  pot += -dr*HubbCi*gconst*MHubbbulge*log(radius/RHubbbulge + r2term);

  fx += -x*dr*Hubbforce;
  fy += -y*dr*Hubbforce;
  fz += -z*dr*Hubbforce;

  pots = - dr*HubbCi2*gconst*MHubbbulge2*log(radius/RHubbbulge2 + r2term2)
         + dr*HubbCi *gconst*MHubbbulge *log(radius/RHubbbulge  + r2term );
  fxs  = - x*dr*Hubbforce2 + x*dr*Hubbforce;
  fys  = - y*dr*Hubbforce2 + y*dr*Hubbforce;
  fzs  = - z*dr*Hubbforce2 + z*dr*Hubbforce;

  if ( ti <= FullBt)
  {

     pot += potb*ti/FullBt + pots*(1.0-ti/FullBt);
     fx  += fxb *ti/FullBt + fxs *(1.0-ti/FullBt);
     fy  += fyb *ti/FullBt + fys *(1.0-ti/FullBt);
     fz  += fzb *ti/FullBt + fzs *(1.0-ti/FullBt);

  } else 

  {

    pot += potb;
    fx  += fxb; 
    fy  += fyb;
    fz  += fzb;

  }

  return;

}
