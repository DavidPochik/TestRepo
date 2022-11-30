//==============================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//==============================================================================
//! \file parker_my_IC.cpp
//  \brief Problem generator for Parker wind.
//


#define COMP_DT 1    // Compute hydro time-steps and save to uov

// C/C++ headers
#include <algorithm>  // min, max
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>   // std::numeric_limits<float>::epsilon()
#include <sstream>
#include <string>
// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"                  // Globals
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../scalars/scalars.hpp"

// Configuration checking

// Static Variables
// critical radius, inner radius, central gravitational parameter,
// isothermal sound speed, initial phi-velocity, initial field strength,
// tilt angle, frame angular velocity
static Real r_0, inv_r2, rho_0, rho_f, v_f, p_f, mu, Ye, Na, p_0, p0, Tc, Th, tau, T_0, dpdd_0, B_0, alpha, T_floor;
static Real Coeff_nu_0, Coeff_nubar_0, t_L_0, t_L_1, t_coeff, dCoeff_nu, dCoeff_nubar;
static const Real float_eps = std::numeric_limits<Real>::epsilon();
static int rows; //number of rows in the IC file
static int IDT1, IDT2, IDT3, IDT4,IDT5,IDT6;
static bool use_IC_file, set_T_at_r0;
Real fact, vfactor;
static int ye_index, IYE;
static Real ye_const, ye_0;
Real* edens;
Real* efrac;
Real* edens0;
Real* Ye0;

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

// Heating/cooling function
Real qdotQW(Real temp, Real x, Real ye, Real time) {
  // temp is in MeV
  // returns heating - cooling in units of MeV s^{-1} g^{-1} units

  // smoothly transition from 0 at t<=t_L_0 to 1 at t>=t_L_1
  Real f = (time <= t_L_0) ? 0.0 :
           ((time >= t_L_1) ? 1.0 : SQR(std::sin((time - t_L_0) * t_coeff)));
  Real Coeff_nu = Coeff_nu_0 + f * dCoeff_nu;
  Real Coeff_nubar = Coeff_nubar_0 + f * dCoeff_nubar;
  // Heating; multiplied 1e12 becasue r is units of 1e6 cm (see Qian and Woosely)
  Real out=1e12*9.65*Na*((1.0-ye)*Coeff_nu + ye*Coeff_nubar)*(1.0-x)*inv_r2;
  out -= 2.27*Na*std::pow(temp,6); // cooling
  if(temp<Tc) {
    out *= std::exp(-0.5/temp); //alpha particles turnoff heating term for T<0.5 MeV
  }
//  if(temp>Th){
//    out *= std::exp(-tau);
//  }
  return out;
}

// Temperature at which heating=cooling
inline Real zeroQ_temp(Real x, Real ye, Real time=0.0) {
  // smoothly transition from 0 at t<=t_L_0 to 1 at t>=t_L_1
  Real f = (time <= t_L_0) ? 0.0 :
           ((time >= t_L_1) ? 1.0 : SQR(std::sin((time - t_L_0) * t_coeff)));
  Real Coeff_nu = Coeff_nu_0 + f * dCoeff_nu;
  Real Coeff_nubar = Coeff_nubar_0 + f * dCoeff_nubar;
  return std::pow(1e12*9.65/2.27*((1.-ye)*Coeff_nu+ye*Coeff_nubar)*(1.-x)*inv_r2,
                  1./6.) / 8.6173e-11;
}


// exists_test1 from https://stackoverflow.com/a/12774387/2275975
inline bool exists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

//Heating/cooling source term
void heat_cool(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
               const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
               AthenaArray<Real> &cons_scalar);

// Boundary Condition
void InflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh);
void InflowMdotInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                       FaceField &b, Real time, Real dt, int is, int ie, int js,
                       int je, int ks, int ke, int ngh);
void HSEInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                FaceField &b, Real time, Real dt, int is, int ie, int js,
                int je, int ks, int ke, int ngh);
void HSE2InnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                 FaceField &b, Real time, Real dt, int is, int ie, int js,
                 int je, int ks, int ke, int ngh);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh);
void InflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh);

//------------------------------------------------------------------------------
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.

void Mesh::InitUserMeshData(ParameterInput *pin) {
  EnrollUserExplicitSourceFunction(heat_cool);
  // set outer BC
  if (pin->GetString("mesh", "ox1_bc").compare("user") == 0) {
    if (Globals::my_rank == 0)
      printf("Using USER outer BC.\n");
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, InflowOuterX1);
  }
  // set inner BC
  int inner_BC_choice = pin->GetOrAddInteger("problem", "inner_BC_choice", 0);
  if (inner_BC_choice == 0) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InflowInnerX1);
  } else if (inner_BC_choice == 1) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, InflowMdotInnerX1);
  } else if (inner_BC_choice == 2) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HSEInnerX1);
  } else if (inner_BC_choice == 3) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, HSE2InnerX1);
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
        << "Invalid inner BC choice (" << inner_BC_choice << ")." << std::endl;
    ATHENA_ERROR(msg);
  }
  printf("Using USER inner BC %d.\n", inner_BC_choice);

  mu = pin->GetReal("problem","GM");
  rho_0 = pin->GetReal("problem","rho_0");
  rho_f = pin->GetReal("problem","rho_f");
  v_f = pin->GetReal("problem","v_f");
  p_f = pin->GetReal("problem","p_f");
  r_0 = pin->GetReal("mesh","x1min");
  p0 = pin->GetReal("problem","p_0");
  inv_r2 = std::pow(r_0, -2);
  Na = pin->GetReal("problem","Na");
  Ye = pin->GetReal("problem","Ye");
  Tc = pin->GetReal("problem","T_cutoff");
  Th = pin->GetReal("problem","T_high");
  tau = pin->GetReal("problem","tau_h");
  B_0 = pin->GetReal("problem","B_0");
  B_0 = B_0/(std::pow(4.0*PI,0.5)); //convert to Lorentz-Heaviside units
  alpha = pin->GetReal("problem","alpha");
  set_T_at_r0 = pin->GetOrAddBoolean("problem", "set_T_at_r0", false);
  // Passive scalar quantities
  ye_index = pin->GetInteger("hydro","eos_ye_index");
  IYE      = pin->GetInteger("hydro","prim_ye_index");
  ye_const = pin->GetReal("hydro","ye_constant");
  ye_0     = pin->GetReal("hydro","Ye_0");
  // final limonosity/energies
  Real L_nu = pin->GetReal("problem","L_nu");
  Real L_nubar = pin->GetReal("problem","L_nubar");
  Real eps_nu = pin->GetReal("problem","eps_nu");
  Real eps_nubar = pin->GetReal("problem","eps_nubar");
  Coeff_nu_0 = L_nu * SQR(eps_nu);
  Coeff_nubar_0 = L_nubar * SQR(eps_nubar);
  // finial limonosity/energies
  L_nu = pin->GetOrAddReal("problem","L_nu_f",L_nu);
  L_nubar = pin->GetOrAddReal("problem","L_nubar_f",L_nubar);
  eps_nu = pin->GetOrAddReal("problem","eps_nu_f",eps_nu);
  eps_nubar = pin->GetOrAddReal("problem","eps_nubar_f",eps_nubar);
  Real coeff = L_nu * SQR(eps_nu);
  dCoeff_nu = coeff - Coeff_nu_0;
  coeff = L_nubar * SQR(eps_nubar);
  dCoeff_nubar = coeff - Coeff_nubar_0;
  Real inf = std::numeric_limits<Real>::infinity();
  t_L_0 = pin->GetOrAddReal("problem","l_transition_start",inf);
  t_L_1 = pin->GetOrAddReal("problem","l_transition_end",inf);
  if (t_L_1 < inf && t_L_1 <= t_L_0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
        << "l_transition_end <= l_transition_start" << std::endl;
    ATHENA_ERROR(msg);
  }
  t_coeff = 0.5 * PI / (t_L_1 - t_L_0);
  T_floor = pin->GetOrAddReal("hydro", "T_floor", float_eps);

  // Parse IC choice
  std::string file;
  bool has_file = pin->DoesParameterExist("problem", "file");
  if (has_file) {
    file = pin->GetString("problem", "file");
  }
  bool use_IC_specified = pin->DoesParameterExist("problem", "use_IC_file");
  use_IC_file = pin->GetOrAddBoolean("problem", "use_IC_file", has_file);

  if (use_IC_specified && use_IC_file) {
    if (!has_file) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
          << "No IC file specified in input file." << std::endl;
      ATHENA_ERROR(msg);
    }
    if (!exists(file)) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Mesh::InitUserMeshData" << std::endl
          << "Specified IC file " << file << "does not exits." << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  if (has_file) {
    if (!exists(file)) {
      use_IC_file = false;
      if (Globals::my_rank == 0) {
        std::cout << "Unable to locate IC file " << file << ", reverting to code IC."
                  << std::endl;
      }
    }
  }

  // Read ICs from data file
  if (use_IC_file) {
    rows = pin->GetInteger("problem", "rows");
    int cols = pin->GetInteger("problem", "cols");
    int col_rho = pin->GetInteger("problem", "col_rho");
    int col_v = pin->GetInteger("problem", "col_v");
    int col_T = pin->GetInteger("problem", "col_T");

    // Prepare arrays to hold profile
    AllocateRealUserMeshDataField(7);
    ruser_mesh_data[0].NewAthenaArray(rows);
    ruser_mesh_data[1].NewAthenaArray(rows);
    ruser_mesh_data[2].NewAthenaArray(rows);
    ruser_mesh_data[3].NewAthenaArray(rows);
    AthenaArray<Real>& r_in{ruser_mesh_data[0]};
    AthenaArray<Real>& rho_in{ruser_mesh_data[1]};
    AthenaArray<Real>& v_in{ruser_mesh_data[2]};
    AthenaArray<Real>& T_in{ruser_mesh_data[3]};

    if (Globals::my_rank == 0)
      std::cout<< "Using IC file: " << file << "\n";
    std::string line;
    std::ifstream stream;

    stream.open(file);
    Real s_vals[cols];

    for (int n = 0; n < rows; ++n) {
      std::getline(stream, line);
      std::string word;
      std::stringstream iss(line);
      int m=0;
      while (iss >> word) {
        s_vals[m]=std::stof(word);
        m=m+1;
      }
      //std::cout<<line;
      r_in(n)=s_vals[0];
      rho_in(n) = s_vals[col_rho+1];
      v_in(n) = s_vals[col_v+1];
      T_in(n) = s_vals[col_T+1];
      if (Globals::my_rank == 0) {
       // std::cout<<r_in(n)<<" ,"<<rho_in(n)<<" , "<<v_in(n)<<" , "<<T_in(n);
       // std::cout<<"\n";
      }
    }

  }
  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  pscalars->s(ye_index,0,0,0) = ye_0*rho_0;
  pscalars->r(ye_index,0,0,0) = pscalars->s(ye_index,0,0,0)/rho_0;
  edens0 = &pscalars->s(ye_index,0,0,0);
  Ye0    = &pscalars->r(ye_index,0,0,0);
  T_0    = zeroQ_temp(0.0, ye_0);
  p_0    = peos->PresFromRhoT(rho_0, T_0, edens0);
  dpdd_0 = peos->AsqFromRhoP(rho_0, p_0, Ye0);

  if (COMP_DT) {
      int i = 0;

      IDT1 = i;
      IDT2 = i + 1;
      IDT3 = i + 2;
      IDT4 = i + 3;
      IDT5 = i + 4;
      IDT6 = i + 5;
      i += 6;

      AllocateUserOutputVariables(i);

      SetUserOutputVariableName(IDT1, "dt1");
      SetUserOutputVariableName(IDT2, "dt2");
      SetUserOutputVariableName(IDT3, "dt3");
      SetUserOutputVariableName(IDT4, "x1flux");
      SetUserOutputVariableName(IDT5, "dflx_vol");
      SetUserOutputVariableName(IDT6, "coord_src1");
  }
}

//------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Parker wind

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // if root processor and zeroth block
 // if ((Globals::my_rank == 0) && (lid == 0)){
 //   std::cout<<"ENTER PGEN NOW NOW, finally right\n";
 // }

  if (use_IC_file) {
    // define references for MeshBlock::ProblemGenerator
    AthenaArray<Real>& r_in{pmy_mesh->ruser_mesh_data[0]};
    AthenaArray<Real>& rho_in{pmy_mesh->ruser_mesh_data[1]};
    AthenaArray<Real>& v_in{pmy_mesh->ruser_mesh_data[2]};
    AthenaArray<Real>& T_in{pmy_mesh->ruser_mesh_data[3]};

    for (int k=ks; k<=ke; k++) {
      // Real phi = pcoord->x3v(k);
      for (int j=js; j<=je; j++) {
        Real theta = pcoord->x2v(j);
     //	std::cout<< "x2v(j) = " << pcoord->x2v(j) << "\n";
        for (int i=is; i<=ie; i++) {
          Real r = pcoord->x1v(i);
          Real r0 = 5e6;
          Real rho, v, temp;

          int index=0;
          Real min=1e15;
          Real diff;

          for (int f=0; f<rows; f++) {
            diff=r-r_in(f);
            if(diff>=0.0) {
              if(diff<min) {
                min=diff;
                index=f;
              }
            }
          }
          //linear interpolation when r values in ICs and simulation are different
          if(r<2.1e6 and rho_0>1.5e12) {
            Real qp= std::pow(r_0/r,20.0);
            rho=rho_0*qp;
            // Real mdot= 4.0*3.14*r*r*rho_in(index)*v_in(index);
            v=(4.34e-4)*2e33/(4.0*3.14*r*r*rho);
          }
          else {
            rho=rho_in(index)+(r-r_in(index))*(rho_in(index+1)-rho_in(index))
                    /(r_in(index+1)-r_in(index));
            v=v_in(index)+(r-r_in(index))*(v_in(index+1)-v_in(index))
                  /(r_in(index+1)-r_in(index));
          }

          temp=T_in(index)+(r-r_in(index))*(T_in(index+1)-T_in(index))
               /(r_in(index+1)-r_in(index));
          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = v * rho;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          pscalars->s(ye_index,k,j,i) = ye_const*rho;
          pscalars->r(ye_index,k,j,i) = pscalars->s(ye_index,k,j,i)/rho;
          edens = &pscalars->s(ye_index,k,j,i);
          efrac = &pscalars->r(ye_index,k,j,i);
          if (NON_BAROTROPIC_EOS) {
            if (GENERAL_EOS) {
              Real pressure = peos->PresFromRhoT(rho, temp, edens);
              phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(rho, pressure, efrac);
            }
            phydro->u(IEN,k,j,i)+=0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       +SQR(phydro->u(IM3,k,j,i)))/rho;
          }
        }
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    // if root processor and zeroth block
    if ((Globals::my_rank == 0) && (lid == 0)){
      std::cout<<"YES ENTER B field\n";
    }
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
    int jl=js; int ju=je+1;
    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary

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
    // Calculate cell-centered magnetic field
    AthenaArray<Real> bb;
    bb.NewAthenaArray(3, ke+1, je+1, ie+NGHOST+1);
    pfield->CalculateCellCenteredField(pfield->b, bb, pcoord, is-NGHOST, ie+NGHOST, js,
                                       je, ks, ke);

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real& bcc1 = bb(IB1,k,j,i);
          Real& bcc2 = bb(IB2,k,j,i);
          Real& bcc3 = bb(IB3,k,j,i);

          phydro->u(IEN,k,j,i) += 0.5*(SQR(bcc1)+SQR(bcc2)+SQR(bcc3));
         }
      }
    }
    a1.DeleteAthenaArray();
    a2.DeleteAthenaArray();
    a3.DeleteAthenaArray();
    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
    bb.DeleteAthenaArray();
  } // end if MAGNETIC_FIELDS_ENABLED
  // if root processor and last block
 // if ((Globals::my_rank == 0) && (lid == pmy_mesh->nblocal - 1)){
 //   std::cout<<"EXIT PGEN NOW NOW, finally right lol\n";
 // }
}


// Source Terms
void heat_cool(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
               const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
               AthenaArray<Real> &cons_scalar) {
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real r = pmb->pcoord->x1v(i);
        Real p = prim(IPR,k,j,i);
        Real rho = prim(IDN,k,j,i);
        Real elecfrac = prim_scalar(ye_index,k,j,i);
        edens = &cons_scalar(ye_index,k,j,i);
        Real temp = pmb->peos->TFromRhoP(rho, p, edens) * 8.6173e-11; //T in MeV
//        std::cout << "ne            = " << cons(IYE,k,j,i) << " cm^-3" << std::endl;
//        std::cout << "Ye            = " << prim(IYE,k,j,i) << std::endl;
//        std::cout << "cons_scalar   = " << cons_scalar(ye_index,k,j,i) << " cm^-3" << std::endl;
//        std::cout << "prim_scalar   = " << prim_scalar(ye_index,k,j,i) << std::endl;
//        std::cout << "----------------------------------------------" << std::endl;
//        std::cout << "rho = " << rho << " g/cm^3" << std::endl;
//        std::cout << "---------------------" << std::endl;
//        std::cout << "Ye (heat_cool) = " << elecfrac << std::endl;

        Real x = std::sqrt(1.0-(r_0*r_0)/(r*r));

        Real qdot = qdotQW(temp, x, elecfrac, time); //MeV s^{-1} g^{-1} units
        Real de = dt*prim(IDN,k,j,i) * qdot * 1.6022e-6; //removed tanh
        cons(IEN,k,j,i) += de;
      }
    }
  }
  return;
}


// Inflow Boundary Condition
void InflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
    // Real phi = pco->x3v(k);
    for (int j=js; j<=je; ++j) {
      Real theta = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        // Real r = pco->x1v(is);
        prim(IDN,k,j,is-i) = rho_0;
        prim(IVX,k,j,is-i) = std::max(prim(IVX,k,j,is), 0.0);
        prim(IVY,k,j,is-i) = 0.0;
        prim(IVZ,k,j,is-i) = 0.0;
        pmb->pscalars->r(ye_index,k,j,is-i) = ye_0;
//	prim(IPR,k,j,is-i) = p0;
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

// Inflow Boundary Condition
void InflowMdotInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
    // Real phi = pco->x3v(k);
    for (int j=js; j<=je; ++j) {
      Real theta = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        prim(IDN,k,j,is-i) = rho_0;
        Real v = prim(IVX,k,j,is) * SQR(pco->x1v(is) / pco->x1v(is-i));
        prim(IVX,k,j,is-i) = std::max(v, 0.0);
        prim(IVY,k,j,is-i) = 0.0;
        prim(IVZ,k,j,is-i) = 0.0;
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

void HSEInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                FaceField &b, Real time, Real dt, int is, int ie, int js,
                int je, int ks, int ke, int ngh) {

  Real r0 = pco->x1v(is);
  for (int k=ks; k<=ke; ++k) {
    // Real phi = pco->x3v(k);
    for (int j=js; j<=je; ++j) {
      Real theta = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        Real r = pco->x1v(is-i);
        prim(IDN,k,j,is-i) = rho_0 * std::exp((r0 - r) * mu / (dpdd_0 * r * r0));
        Real v = prim(IVX,k,j,is) * prim(IDN,k,j,is) / prim(IDN,k,j,is-i) * SQR(r0 / r);
        prim(IVX,k,j,is-i) = std::max(v, 0.0);
        prim(IVY,k,j,is-i) = 0.0;
        prim(IVZ,k,j,is-i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          prim(IPR,k,j,is-i) = pmb->peos->PresFromRhoT(prim(IDN,k,j,is-i), T_0);
        }
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

void HSE2InnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                 FaceField &b, Real time, Real dt, int is, int ie, int js,
                int je, int ks, int ke, int ngh) {
  AthenaArray<Real> out1;
  out1.NewAthenaArray(7);
  Real r0 = pco->x1v(is);
  for (int k=ks; k<=ke; ++k) {
    // Real phi = pco->x3v(k);
    for (int j=js; j<=je; ++j) {
      Real theta = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        Real r = pco->x1v(is-i);
        prim(IDN,k,j,is-i) = rho_0 * std::exp((r0 - r) * mu / (dpdd_0 * r * r0));
        prim(IVX,k,j,is-i) = std::max(prim(IVX,k,j,is), 0.0);
        prim(IVY,k,j,is-i) = 0.0;
        prim(IVZ,k,j,is-i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          prim(IPR,k,j,is-i) = pmb->peos->PresFromRhoT(prim(IDN,k,j,is-i), T_0);
        }
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
      for (int i=1; i<=ngh; ++i) {
        Real rgh = pco->x1v(ie+i);
        Real ratio = re / rgh;
        prim(IDN,k,j,ie+i) = prim(IDN,k,j,ie) * SQR(ratio);
        prim(IVX,k,j,ie+i) = std::max(prim(IVX,k,j,ie), 0.0);
        prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie) * ratio;
        prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie) * ratio;
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

// Inflow Boundary Condition at outer boundary
void InflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt, int is, int ie, int js,
                   int je, int ks, int ke, int ngh) {
 // Real re = pco->x1v(ie);
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real theta = pco->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        prim(IDN,k,j,ie+i) = rho_f;
        prim(IVX,k,j,ie+i) = v_f;
        prim(IPR,k,j,ie+i) = p_f;
        prim(IVY,k,j,ie+i) = 0.0;
        prim(IVZ,k,j,ie+i) = 0.0;
        pmb->pscalars->r(ye_index,k,j,ie+i)=ye_const;
   //     if (NON_BAROTROPIC_EOS)
   //       prim(IPR,k,j,ie+i) = prim(IPR,k,j,ie); // Constant Pressure
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
  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {

      for(int i=is; i<=ie; i++) {
//        pscalars->s(ye_index,k,j,i)=ye_const*phydro->u(IDN,k,j,i);
//        pscalars->r(ye_index,k,j,i)=pscalars->s(ye_index,k,j,i)/phydro->u(IDN,k,j,i);
        edens = &pscalars->s(ye_index,k,j,i);
        efrac = &pscalars->r(ye_index,k,j,i);
//        std::cout << "edens = " << edens << std::endl;
//        std::cout << "efrac (UserWorkBeforeOutput) = " << phydro->u(IYE,k,j,i)/phydro->u(IDN,k,j,i) << std::endl;
        Real temp=peos->TFromRhoP(phydro->w(IDN,k,j,i),phydro->w(IPR,k,j,i),edens);
        // Real egas= peos->EgasFromRhoP(phydro->w(IDN,k,j,i),phydro->w(IPR,k,j,i));
        // Real temp2= peos->TFromRhoEgas(phydro->w(IDN,k,j,i), egas);

        user_out_var(0,k,j,i) = temp;
        Real r = pcoord->x1v(i);

        temp *= 8.6173e-11; //convert to MeV
        Real x=std::pow((1.0-(r_0*r_0)/(r*r)),0.5);
        Real qdot = qdotQW(temp, x, pscalars->r(ye_index,k,j,i), pmy_mesh->time) * 1.6022e-6; // in ergs/s/g
        user_out_var(1,k,j,i) = qdot;
        user_out_var(2,k,j,i) = std::sqrt(peos->AsqFromRhoP(phydro->w(IDN,k,j,i),
                                                            phydro->w(IPR,k,j,i),efrac));
      }
    }
  }
}
