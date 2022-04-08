//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file parker_wind.cpp
//  \brief Problem generator for Parker wind.
//

#define COMP_DIV_B 1
#define COMP_DT 1    // Compute hydro time-steps and save to uov

// C/C++ headers
#include <algorithm>  // min, max
#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"

// Configuration checking

// Static Variables
// critical radius, inner radius, central gravitational parameter,
// isothermal sound speed, initial phi-velocity, initial field strength,
// tilt angle, frame angular velocity
static Real r_crit, r_0, rho_0, mu, a_iso, vr_0, vphi_0, B_0, alpha, Omega_0, vpow,
       Pconst, p_0;
static int IDIVB, IDT1, IDT2, IDT3;

// User-Defined History Output
Real MdotHst(MeshBlock *pmb, int iout);
Real JdotHst(MeshBlock *pmb, int iout);
Real EdotHst(MeshBlock *pmb, int iout);

// Rotarional Source Terms
void RotFrame(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
              AthenaArray<Real> &cons);

// Boundary Condition
void InflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh);

// Profile Calculations
static Real HydrostaticDensityProfile(const Real r);
static void HSEprof(const Real r, Real &rho, Real &v, Real &p);
static Real ConservedRotationProfile(const Real r, const Real theta);

// Vector Potential
static Real A3(const Real x1, const Real x2, const Real x3);
static Real A2(const Real x1, const Real x2, const Real x3);
static Real A1(const Real x1, const Real x2, const Real x3);

//------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.

void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real rot, mag, vesc_0;
  EnrollUserExplicitSourceFunction(RotFrame);
  if (pin->GetString("mesh", "ox1_bc").compare("user") == 0) {
      EnrollUserBoundaryFunction(BoundaryFace::outer_x1, OutflowOuterX1);
  }
  EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InflowInnerX1);
  AllocateUserHistoryOutput(4);
  EnrollUserHistoryOutput(0,MdotHst,"Mdot");
  EnrollUserHistoryOutput(1,JdotHst,"Jdot");
  EnrollUserHistoryOutput(2,EdotHst,"Edot");
  mu = pin->GetReal("problem","GM");
  rho_0 = pin->GetReal("problem","rho_0");
  r_0 = pin->GetReal("mesh","x1min");
  a_iso = pin->GetReal("hydro","iso_sound_speed");
  p_0 = pin->GetOrAddReal("problem", "p_0", SQR(a_iso) * rho_0);
  vpow = pin->GetOrAddReal("problem","v_power", 0.5);
  rot = pin->GetOrAddReal("problem","rot",0.0);
  mag = pin->GetOrAddReal("problem","mag",0.0);
  alpha = pin->GetOrAddReal("problem","alpha",0.0);
  vesc_0 = std::sqrt(2*mu / r_0);
  B_0 = std::sqrt(rho_0) * mag * vesc_0;
  vphi_0 = rot * vesc_0;
  Omega_0 = vphi_0 / r_0;
  r_crit = mu / (2*a_iso*a_iso);
  vr_0 = a_iso * std::pow(r_0/r_crit, vpow);
  Pconst = rho_0 * mu/(r_0 * (1.0 + vpow)) + SQR(a_iso) * rho_0 + p_0;
  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  if ((COMP_DIV_B & MAGNETIC_FIELDS_ENABLED) || COMP_DT) {
    int i = 0;
    if (COMP_DIV_B & MAGNETIC_FIELDS_ENABLED) {
      IDIVB = i;
      i ++;
    }
    if (COMP_DT) {
      IDT1 = i;
      IDT2 = i + 1;
      IDT3 = i + 2;
      i += 3;
    }
    AllocateUserOutputVariables(i);
    if (COMP_DIV_B & MAGNETIC_FIELDS_ENABLED) SetUserOutputVariableName(IDIVB, "divB");
    if (COMP_DT) {
      SetUserOutputVariableName(IDT1, "dt1");
      SetUserOutputVariableName(IDT2, "dt2");
      SetUserOutputVariableName(IDT3, "dt3");
    }
  }
}

//------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Parker wind

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Read problem parameters
  Real r,theta,phi,gm1;
  if (NON_BAROTROPIC_EOS && !GENERAL_EOS)
    gm1 = peos->GetGamma() - 1.0;
  for (int k=ks; k<=ke; k++) {
    phi = pcoord->x3v(k);
    for (int j=js; j<=je; j++) {
      theta = pcoord->x2v(j);
      for (int i=is; i<=ie; i++) {
        r = pcoord->x1v(i);
        Real rho, v, P;
        HSEprof(r, rho, v, P);
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = v * rho;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = ConservedRotationProfile(r,theta) * rho;
        if (NON_BAROTROPIC_EOS) {
          if (GENERAL_EOS) {
            phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(rho, P);
          } else {
            phydro->u(IEN,k,j,i) = P / gm1;
          }
          phydro->u(IEN,k,j,i) += (SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                   +SQR(phydro->u(IM3,k,j,i)))/rho;
        }
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    AthenaArray<Real> a1,a2,a3;
    int nx1 = (ie-is)+1 + 2*(NGHOST);
    int nx2 = (je-js)+1 + 2*(NGHOST);
    int nx3 = (ke-ks)+1 + 2*(NGHOST);
    a1.NewAthenaArray(nx3,nx2,nx1);
    a2.NewAthenaArray(nx3,nx2,nx1);
    a3.NewAthenaArray(nx3,nx2,nx1);

    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
          a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
          a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
        }
      }
    }

    // Initialize interface fields
    AthenaArray<Real> area,len,len_p1;
    area.NewAthenaArray(nx1);
    len.NewAthenaArray(nx1);
    len_p1.NewAthenaArray(nx1);

    // for 1,2,3-D
    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      int jl=js; int ju=je+1;
      if ((pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) ||
        (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge))
        jl=js+1;
      if ((pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) ||
        (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
        ju=je;
      for (int j=jl; j<=ju; ++j) {
        pcoord->Face2Area(k,j,is,ie,area);
        pcoord->Edge3Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = -1.0*(len(i+1)*a3(k,j,i+1) - len(i)*a3(k,j,i))/area(i);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        pcoord->Face3Area(k,j,is,ie,area);
        pcoord->Edge2Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = (len(i+1)*a2(k,j,i+1) - len(i)*a2(k,j,i))/area(i);
        }
      }
    }

    // for 2D and 3D
    if (block_size.nx2 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge3Length(k,j  ,is,ie+1,len);
          pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = (len_p1(i)*a3(k,j+1,i) - len(i)*a3(k,j,i))/area(i);
          }
        }
      }
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face3Area(k,j,is,ie,area);
          pcoord->Edge1Length(k,j  ,is,ie,len);
          pcoord->Edge1Length(k,j+1,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) -= (len_p1(i)*a1(k,j+1,i) - len(i)*a1(k,j,i))/area(i);
          }
        }
      }
    }
    // for 3D only
    if (block_size.nx3 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge2Length(k  ,j,is,ie+1,len);
          pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) -= (len_p1(i)*a2(k+1,j,i) - len(i)*a2(k,j,i))/area(i);
          }
        }
      }
      for (int k=ks; k<=ke; ++k) {
        // reset loop limits for polar boundary
        int jl=js; int ju=je+1;
        if ((pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) ||
          (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge))
          jl=js+1;
        if ((pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) ||
          (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge))
          ju=je;
        for (int j=jl; j<=ju; ++j) {
          pcoord->Face2Area(k,j,is,ie,area);
          pcoord->Edge1Length(k  ,j,is,ie,len);
          pcoord->Edge1Length(k+1,j,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) += (len_p1(i)*a1(k+1,j,i) - len(i)*a1(k,j,i))/area(i);
          }
        }
      }
    }

    a1.DeleteAthenaArray();
    a2.DeleteAthenaArray();
    a3.DeleteAthenaArray();
    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
}
}

// Computes Hydrostatic Density profile
static Real HydrostaticDensityProfile(const Real r) {
  Real grav, mach;
  Real ir, ir_0;
  ir = 1.0/r; // Inverse r
  ir_0 = 1.0/r_0; // Inverse r_0
  mach = (vphi_0*vphi_0)/(a_iso*a_iso); // Rotational Mach number squared
  grav = std::exp((mu/(a_iso*a_iso)) * (ir - ir_0)); // Gravitational HSE Term
  return rho_0 * grav;
}

// Computes Hydrostatic profile
static void HSEprof(const Real r, Real &rho, Real &v, Real &p) {
  Real rrat = r / r_0;
  v = vr_0 * std::pow(rrat, vpow);
  if (NON_BAROTROPIC_EOS) {
    Real rho_pow = -vpow - 2.0;
    rho = rho_0 * std::pow(rrat, rho_pow);
    p = rho / (1.0 + rho_pow) * mu / r - rho * SQR(v) + Pconst / SQR(rrat);
  } else {
    Real grav;
    Real ir, ir_0;
    ir = 1.0/r; // Inverse r
    ir_0 = 1.0/r_0; // Inverse r_0
    grav = std::exp((mu/(a_iso*a_iso)) * (ir - ir_0)); // Gravitational HSE Term
    rho = rho_0 * grav;
  }
}

static Real ConservedRotationProfile(const Real x1, const Real x2) {
  return Omega_0 * std::sin(x2) * x1 * (std::pow((r_0/x1),2) - 1);
}

// Vector Potential
static Real A3(const Real x1, const Real x2, const Real x3) {
  Real a3 = 0.5 * B_0 * r_0 * std::pow(r_0/x1,2) *
    (std::sin(x2)*std::cos(alpha) - std::cos(x2)*std::cos(x3)*std::sin(alpha));
  return a3;
}

static Real A2(const Real x1, const Real x2, const Real x3) {
  Real a2 = -0.5 * B_0 * r_0 * std::pow(r_0/x1,2) * std::sin(x3) * std::sin(alpha);
  return a2;
}

static Real A1(const Real x1, const Real x2, const Real x3) {
  Real a1=0.0;
  return a1;
}

// User-Defined History Output
Real MdotHst(MeshBlock *pmb, int iout) {
  AthenaArray<Real> vol;
  vol.NewAthenaArray(pmb->ncells1);
  Real mdot = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real vel1 = pmb->phydro->w(IM1,k,j,i);
        Real rho = pmb->phydro->w(IDN,k,j,i);
        mdot += vel1*rho*vol(i);
      }
    }
  }
  return mdot;
}

Real JdotHst(MeshBlock *pmb, int iout) {
  AthenaArray<Real> vol;
  vol.NewAthenaArray(pmb->ncells1);
  Real jdot = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      Real sinx2 = std::sin(pmb->pcoord->x2v(j));
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real rho = pmb->phydro->w(IDN,k,j,i);
        Real vel1 = pmb->phydro->w(IM1,k,j,i);
        Real vel3 = pmb->phydro->w(IM3,k,j,i) + Omega_0*x1*sinx2;
        jdot += x1 * vel1 * vel3 * rho * vol(i);
        if (MAGNETIC_FIELDS_ENABLED) {
          Real Bcc1 = 0.5 * (pmb->pfield->b.x1f(k,j,i) + pmb->pfield->b.x1f(k,j,i+1));
          Real Bcc3 = 0.5 * (pmb->pfield->b.x3f(k,j,i) + pmb->pfield->b.x3f(k+1,j,i));
          jdot -= x1 * Bcc1 * Bcc3 * vol(i);
        }
      }
    }
  }
  return jdot;
}

Real EdotHst(MeshBlock *pmb, int iout) {
  AthenaArray<Real> vol;
  vol.NewAthenaArray(pmb->ncells1);
  Real edot = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      Real sinx2 = std::sin(pmb->pcoord->x2v(j));
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, vol);
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real rho = pmb->phydro->w(IDN,k,j,i);
        Real vel1 = pmb->phydro->w(IM1,k,j,i);
        Real vel2 = pmb->phydro->w(IM2,k,j,i);
        Real vel3 = pmb->phydro->w(IM3,k,j,i) + Omega_0*x1*sinx2;
        edot += rho * std::pow(vel1,3) * vol(i); // Kinetic Flux
        if (MAGNETIC_FIELDS_ENABLED) {
          Real Bcc1 = 0.5 * (pmb->pfield->b.x1f(k,j,i) + pmb->pfield->b.x1f(k,j,i+1));
          Real Bcc2 = 0.5 * (pmb->pfield->b.x2f(k,j,i) + pmb->pfield->b.x2f(k,j+1,i));
          Real Bcc3 = 0.5 * (pmb->pfield->b.x3f(k,j,i) + pmb->pfield->b.x3f(k+1,j,i));
          // Poynting Flux
          edot += (Bcc2*(Bcc2*vel1 - Bcc1*vel2) + Bcc3*(Bcc3*vel1 - Bcc1*vel3)) * vol(i);
        }
      }
    }
  }
  return edot;
}

// Source Terms
void RotFrame(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
              AthenaArray<Real> &cons) {
  // a_cor = -2 * Omega x v
  Real x1,x2,x3;
  Real f = dt * Omega_0;
  Real cen_dt = f * Omega_0;
  Real ff = f*f;
  AthenaArray<Real> &x1flux=pmb->phydro->flux[X1DIR];
  AthenaArray<Real> &x2flux=pmb->phydro->flux[X2DIR];
  // AthenaArray<Real> &x3flux=pmb->phydro->flux[X3DIR];
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    x3 = pmb->pcoord->x3v(k);
    Real dx3 = pmb->pcoord->dx3f(k);
    for (int j=pmb->js; j<=pmb->je; ++j) {
      x2 = pmb->pcoord->x2v(j);
      Real sinx2 = std::sin(x2);
      Real cosx2 = std::cos(x2);
      Real sin2x2 = std::sin(2 * x2);
      Real cos2x2 = std::cos(2 * x2);
      Real sinx2l = std::sin(pmb->pcoord->x2f(j));
      Real sinx2r = std::sin(pmb->pcoord->x2f(j+1));
      Real cosx2l = std::cos(pmb->pcoord->x2f(j));
      Real cosx2r = std::cos(pmb->pcoord->x2f(j+1));
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        x1 = pmb->pcoord->x1v(i);
        Real vol = pmb->pcoord->GetCellVolume(k,j,i);

        // Coriolis
        Real mom1 = cons(IM1,k,j,i);
        Real mom2 = cons(IM2,k,j,i);
        Real mom3 = cons(IM3,k,j,i);
        cons(IM1,k,j,i) += ff*(mom1*cos2x2 - mom2*sin2x2) + 2*f*mom3*sinx2;
        cons(IM1,k,j,i) /= (1+ff);
        cons(IM2,k,j,i) += -ff*(mom1*sin2x2 + mom2*cos2x2) + 2*f*mom3*cosx2;
        cons(IM2,k,j,i) /= (1+ff);
        cons(IM3,k,j,i) -= ff*mom3 + 2*f*(mom1*sinx2 + mom2*cosx2);
        cons(IM3,k,j,i) /= (1+ff);

        // Centrifugal
        Real dr4 = std::pow(pmb->pcoord->x1f(i+1),4) - std::pow(pmb->pcoord->x1f(i),4);
        Real cen_0 = 1./12. * dx3 / vol * cen_dt * dr4;
        Real src1 = .25*cen_0*(-9*(cosx2r-cosx2l)
                               -3*(cosx2r*SQR(sinx2r)-cosx2l*SQR(sinx2l))
                               +std::pow(cosx2r,3) -std::pow(cosx2l,3));
        cons(IM1,k,j,i) += src1 * prim(IDN,k,j,i);
        Real src2 = cen_0*(std::pow(sinx2r,3)-std::pow(sinx2l,3));
        cons(IM2,k,j,i) += src2 * prim(IDN,k,j,i);

        // Energy
        if (NON_BAROTROPIC_EOS) {
          // dE = a . v * dt
          cons(IEN,k,j,i) += src1 * 0.5*(x1flux(IDN,k,j,i) + x1flux(IDN,k,j,i+1));
          if (pmb->block_size.nx2 > 1) {
            cons(IEN,k,j,i) += src2 * 0.5*(x2flux(IDN,k,j,i) + x2flux(IDN,k,j+1,i));
          } else {
            cons(IEN,k,j,i) += src2 * prim(IDN,k,j,i) * prim(IVX,k,j,i);
          }
          // (src3=0)     += src3 * 0.5*(x3flux(IDN,k,j,i) + x3flux(IDN,k+1,j,i));
        }
      }
    }
  }
  return;
}

// Inflow Boundary Condition
void InflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh) {
  Real r,theta,phi;
  for (int k=ks; k<=ke; ++k) {
    phi = pco->x3v(k);
    for (int j=js; j<=je; ++j) {
      theta = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        r = pco->x1v(is);
        Real rho = HydrostaticDensityProfile(r);
        prim(IDN,k,j,is-i) = rho;
        prim(IM1,k,j,is-i) = std::max(prim(IM1,k,j,is), 0.0);
        prim(IM2,k,j,is-i) = 0.0;
        prim(IM3,k,j,is-i) = 0.0;
        if (NON_BAROTROPIC_EOS)
          prim(IPR,k,j,is-i) = p_0;
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
        }
      }
    }
  }
}

// Outflow Boundary Condition
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh) {
  Real re = pco->x1v(ie);
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real sinx2 = std::sin(pco->x2v(j));
      for (int i=1; i<=ngh; ++i) {
        Real rgh = pco->x1v(ie+i);
        Real ratio = re / rgh;
        prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie) * SQR(ratio);
        prim(IM1,k,j,ie+i) = std::max(prim(IM1,k,j,ie), 0.0);
        prim(IM2,k,j,ie+i) = prim(IM2,k,j,ie) * ratio;
        prim(IM3,k,j,ie+i) = prim(IM3,k,j,ie)*ratio + Omega_0*rgh*(SQR(ratio) - 1)*sinx2;
        if (NON_BAROTROPIC_EOS)
          prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie); // Constant Pressure
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
        }
      }
    }
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  if (COMP_DIV_B & MAGNETIC_FIELDS_ENABLED) {
    AthenaArray<Real> x1area, x2area, x2area_p1, x3area, x3area_p1, vol, dflx;
    int ncells1 = block_size.nx1+2*NGHOST;
    x1area.NewAthenaArray(ncells1+1);
    x2area.NewAthenaArray(ncells1);
    x2area_p1.NewAthenaArray(ncells1);
    x3area.NewAthenaArray(ncells1);
    x3area_p1.NewAthenaArray(ncells1);
    vol.NewAthenaArray(ncells1);
    dflx.NewAthenaArray(ncells1);
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        // calculate x1-flux divergence
        pcoord->Face1Area(k,j,is,ie+1,x1area);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          dflx(i) = (x1area(i+1)*pfield->b.x1f(k,j,i+1) - x1area(i)*pfield->b.x1f(k,j,i));
        }

        // calculate x2-flux divergence
        if (block_size.nx2 > 1) {
          pcoord->Face2Area(k,j  ,is,ie,x2area   );
          pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(i) += (x2area_p1(i)*pfield->b.x2f(k,j+1,i)
                       -x2area(i)*pfield->b.x2f(k,j,i));
          }
        }

        // calculate x3-flux divergence
        if (block_size.nx3 > 1) {
          pcoord->Face3Area(k  ,j,is,ie,x3area   );
          pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(i) += (x3area_p1(i)*pfield->b.x3f(k+1,j,i)
                       -x3area(i)*pfield->b.x3f(k,j,i));
          }
        }

        // update conserved variables
        pcoord->CellVolume(k,j,is,ie,vol);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          user_out_var(0,k,j,i) = dflx(i)/vol(i);
        }
      }
    }
    x1area.DeleteAthenaArray();
    x2area.DeleteAthenaArray();
    x2area_p1.DeleteAthenaArray();
    x3area.DeleteAthenaArray();
    x3area_p1.DeleteAthenaArray();
    vol.DeleteAthenaArray();
    dflx.DeleteAthenaArray();
  }
  if (COMP_DT) {
    AthenaArray<Real>& w = phydro->w;
    AthenaArray<Real>* bcc;
    AthenaArray<Real>* b_x1f;
    AthenaArray<Real>* b_x2f;
    AthenaArray<Real>* b_x3f;
    if (MAGNETIC_FIELDS_ENABLED) {
      bcc = &pfield->bcc;
      b_x1f = &pfield->b.x1f;
      b_x2f = &pfield->b.x2f;
      b_x3f = &pfield->b.x3f;
    }

    AthenaArray<Real> dt1, dt2, dt3;
    int ncells1 = block_size.nx1+2*NGHOST;
    dt1.NewAthenaArray(ncells1);
    dt2.NewAthenaArray(ncells1);
    dt3.NewAthenaArray(ncells1);
    Real wi[(NWAVE)];

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pcoord->CenterWidth1(k,j,is,ie,dt1);
        pcoord->CenterWidth2(k,j,is,ie,dt2);
        pcoord->CenterWidth3(k,j,is,ie,dt3);
#pragma ivdep
        for (int i=is; i<=ie; ++i) {
          wi[IDN]=w(IDN,k,j,i);
          wi[IVX]=w(IVX,k,j,i);
          wi[IVY]=w(IVY,k,j,i);
          wi[IVZ]=w(IVZ,k,j,i);
          if (NON_BAROTROPIC_EOS)
            wi[IPR]=w(IPR,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED) {
            Real bx = (*bcc)(IB1,k,j,i) + fabs((*b_x1f)(k,j,i)-(*bcc)(IB1,k,j,i));
            wi[IBY] = (*bcc)(IB2,k,j,i);
            wi[IBZ] = (*bcc)(IB3,k,j,i);
            Real cf = peos->FastMagnetosonicSpeed(wi,bx);
            dt1(i) /= (fabs(wi[IVX]) + cf);

            wi[IBY] = (*bcc)(IB3,k,j,i);
            wi[IBZ] = (*bcc)(IB1,k,j,i);
            bx = (*bcc)(IB2,k,j,i) + fabs((*b_x2f)(k,j,i)-(*bcc)(IB2,k,j,i));
            cf = peos->FastMagnetosonicSpeed(wi,bx);
            dt2(i) /= (fabs(wi[IVY]) + cf);

            wi[IBY] = (*bcc)(IB1,k,j,i);
            wi[IBZ] = (*bcc)(IB2,k,j,i);
            bx = (*bcc)(IB3,k,j,i) + fabs((*b_x3f)(k,j,i)-(*bcc)(IB3,k,j,i));
            cf = peos->FastMagnetosonicSpeed(wi,bx);
            dt3(i) /= (fabs(wi[IVZ]) + cf);
          } else {
            Real cs = peos->SoundSpeed(wi);
            dt1(i) /= (fabs(wi[IVX]) + cs);
            dt2(i) /= (fabs(wi[IVY]) + cs);
            dt3(i) /= (fabs(wi[IVZ]) + cs);
          }
          user_out_var(IDT1,k,j,i) = dt1(i) * pmy_mesh->cfl_number;
          user_out_var(IDT2,k,j,i) = dt2(i) * pmy_mesh->cfl_number;
          user_out_var(IDT3,k,j,i) = dt3(i) * pmy_mesh->cfl_number;
        }
      }
    }
    dt1.DeleteAthenaArray();
    dt2.DeleteAthenaArray();
    dt3.DeleteAthenaArray();
  }
}
