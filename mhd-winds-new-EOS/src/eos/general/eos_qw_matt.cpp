//#define MYDEBUG
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file eos_qw.cpp
//! \brief implements functions in class EquationOfState for simple hydrogen EOS
//======================================================================================

// C headers

// C++ headers
#include <algorithm>
#include <cmath>   // sqrt()
#include <fstream>
#include <iostream> // ifstream
#include <limits>   // std::numeric_limits<float>::epsilon()
#include <sstream>
#include <stdexcept> // std::invalid_argument
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../field/field.hpp"
#include "../../parameter_input.hpp"
#include "../eos.hpp"

// Same effect as adding static to everything inside
namespace {
  const Real float_eps = std::numeric_limits<Real>::epsilon();
  // const Real float_1pe = 1.0 + float_eps;
  const Real float_max = std::numeric_limits<Real>::max();
  Real prec = 1e-8;
  Real T_floor, T_ceil, LastTemp;
  Real Ye;
  bool use_T_floor;
  int nmax;
  AthenaArray<Real> EosData;

  const Real third = 1.0 / 3.0;
  const Real c = 3e10;
  const Real k = 1.380649e-16;
  const Real mn = 1.6726e-24;
  const Real hbar = 1.0546e-27;
  const Real c3 = std::pow(k/(hbar*c),3);
  const Real con3 = (11.0*PI*PI/180.0)*c3*k;
}

// EOS data indicies
namespace EOS {
  enum EosIndex {iE=0, idEdT=1, iP=2, idPdT=3, iAsq=4, iT=5, N=6};
}

#ifdef MYDEBUG
namespace compare {
  //returns \eta/pi in terms of Ye. See eqn 6 of QW
  Real eta_by_pi(Real rho, Real T){

    Real a= std::pow(k*T/(hbar*c),3)*(mn/rho)*(PI/3.0);
    Real term= std::pow(9.0*a*a*Ye+1.7321*std::pow((4.0*std::pow(a,6.0)+27.0*std::pow(a,4.0)*Ye*Ye),0.5),1.0/3.0);
    return 0.38157*term/a - 0.87358*a/term;
  }

  //returns (d \eta/pi)/dT in terms of Ye. See eqn 6 of QW
  Real der_etabypi(Real rho, Real T){

    Real a= std::pow(k*T/(hbar*c),3)*(mn/rho)*(PI/3.0);
    Real da_dt= 3.0*T*T*std::pow(k/(hbar*c),3)*(mn/rho)*(PI/3.0);
    Real term= std::pow(9.0*a*a*Ye+1.7321*std::pow((4.0*std::pow(a,6.0)+27.0*std::pow(a,4.0)*Ye*Ye),0.5),1.0/3.0);
    Real dterm_dt= (1.0/3.0)*std::pow(9.0*a*a*Ye+1.7321*std::pow((4.0*std::pow(a,6.0)+27.0*std::pow(a,4.0)*Ye*Ye),0.5),-2.0/3.0)*(18.0*a*Ye*da_dt + 1.7321*0.5*std::pow((4.0*std::pow(a,6.0)+27.0*std::pow(a,4.0)*Ye*Ye),-0.5)*(24.0*std::pow(a,5.0)*da_dt+27.0*4.0*a*a*a*Ye*Ye*da_dt));
    return 0.38157*(dterm_dt/a - term*da_dt/(a*a)) - 0.87358*(da_dt/term -a*dterm_dt/(term*term));


  }
  //! \brief compute gas pressure, see eqn 4 of QW
  Real P_of_rho_T(Real rho, Real T) {
    Real eta = eta_by_pi(rho,T);
    Real con= (11.0*PI*PI/180.0)*std::pow(k/(hbar*c),3)*k*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0);
    return con*std::pow(T,4.0)+ rho*k*T/mn;
  }

  //----------------------------------------------------------------------------------------
  //! \brief compute internal energy density, see eqn 5 of QW
  Real e_of_rho_T(Real rho, Real T) {
    Real eta = eta_by_pi(rho,T);
    Real constant= (11.0*PI*PI/60.0)*std::pow(k/(hbar*c),3)*k*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0);
    Real sp_energy= constant*std::pow(T,4.0) + 1.5*rho*k*T/mn; //energy/g
    return sp_energy;
  }

  //temp derivative of the function
  Real func_derT(Real rho,Real P,Real T) {
    Real eta = eta_by_pi(rho,T);

    Real der_eta= der_etabypi(rho,T);
    Real con3= (11.0*PI*PI/180.0)*std::pow(k/(hbar*c),3)*k;
    Real val2= 4.0*con3*std::pow(T,3.0)*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0)+con3*std::pow(T,4.0)*(30.0*2.0*eta*der_eta/11.0 + 15.0*4.0*std::pow(eta,3)*der_eta/11.0)+ rho*k/mn;
    return val2;
  }

  //----------------------------------------------------------------------------------------
  //! \brief compute adiabatic sound speed squared
  Real asq(Real rho, Real T, Real P) {
    Real eta = eta_by_pi(rho,T);
    Real ctsq=k*T/mn;
    Real der_eta= der_etabypi(rho,T);
    // Real con1= (11.0*PI*PI/60.0)*std::pow(k/(hbar*c),3)*k*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0);
    Real D=(T/rho)*func_derT(rho,P,T);
    Real con3= (11.0*PI*PI/60.0)*std::pow(k/(hbar*c),3)*k;
    Real cv= (1.0/rho)*(4.0*con3*std::pow(T,3.0)*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0)+con3*std::pow(T,4.0)*(30.0*2.0*eta*der_eta/11.0 + 15.0*4.0*std::pow(eta,3)*der_eta/11.0))+ 1.5*k/mn;
    return ctsq+D*D/(cv*T);
  }
}
#endif
//----------------------------------------------------------------------------------------



void QWData(Real rho, Real T, AthenaArray<Real> &OutData){
  Real vol = mn/rho;
  Real T3 = T*T*T;
  Real T4 = T*T3;
  Real a = c3*std::pow(T,3)*vol*(PI/3.0);
  Real a2 = SQR(a);
  Real a4 = SQR(a2);
  Real a6 = a2 * a4;
  Real y2 = SQR(Ye);
  Real b = std::sqrt(4.0*a6+27.0*a4*y2);
  Real term = std::pow(9.0*a2*Ye+1.7321*b, third);
  Real eta = 0.38157*term/a - 0.87358*a/term; // actually eta/pi

  Real eta2 = SQR(eta);
  Real eta3 = eta*eta2;
  Real eta4 = SQR(eta2);

  Real da_dt = 3.0*a/T;
  Real dterm_dt= third*std::pow(term, -2)*(18.0*a*Ye*da_dt
                + 1.7321*0.5*a*(24.0*a4*da_dt+27.0*4.0*a2*y2*da_dt)/b);
  Real der_eta = 0.38157*(dterm_dt/a - term*da_dt/(a*a))
                 - 0.87358*(da_dt/term -a*dterm_dt/(term*term));

  Real con= con3*(1.0+30.0*eta2/11.0 + 15.0*eta4/11.0);
  Real p0 = rho*k*T/mn;
  Real P = con*T4 + p0;
  Real e = 3*con*T4 + 1.5*p0;

  Real der_P= 4.0*T3*con + con3*T4*(30.0*2.0*eta*der_eta/11.0
                                    + 15.0*4.0*eta3*der_eta/11.0);
  Real der_e= 3*der_P+ 1.5*rho*k/mn;
  der_P += rho*k/mn;

  Real ctsq=k*T/mn;
  Real D=(T/rho)*der_P;
  Real cv= der_e / rho;
  Real asq = ctsq+D*D/(cv*T);

  OutData(0) = e;
  OutData(1) = der_e;
  OutData(2) = P;
  OutData(3) = der_P;
  OutData(4) = asq;
  OutData(5) = T;
}

// index = 0 for internal energy; index = 2 for pressure; var = int energy or pressure
void TempInvert(Real rho, Real GuessTemp, Real var, const int index,
                AthenaArray<Real> &OutData) {
  Real T_high;
  if (index == EOS::iP) {
    T_high = std::min(std::pow(var / con3, .25), var * mn / (k * rho));
  } else if (index == EOS::iE) {
    T_high = std::min(std::pow(var / (3.0 * con3), .25), var * mn / (1.5 * k * rho));
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState inversion (TempInvert)"
        << std::endl << "Invalid inversion variable." << std::endl;
    ATHENA_ERROR(msg);
  }
  Real BrakT[] = {T_floor, std::min(T_ceil, T_high)};
  Real BrakVal[] = {0, 0};
  Real InvVar = 1.0 / var;
  Real error = float_max;
  int nlim = nmax;
  Real LastTemp = BrakT[0];
  QWData(rho, LastTemp, OutData);
  BrakVal[0] = OutData(index) * InvVar - 1.0;
  Real LastErr = BrakVal[0];
  QWData(rho, BrakT[1], OutData);
  BrakVal[1] = OutData(index) * InvVar - 1.0;
  Real a = BrakVal[0];
  Real b = BrakVal[1];
  Real delta;
#ifdef MYDEBUG1
  printf("%d: %.16e, %.16e\n", index, var, rho);
  int mode = 0;
#endif
  while (std::abs(error) > prec) {
    if (BrakVal[0] > 0) {//}* BrakVal[1] > 0) {
      QWData(rho, BrakT[0], OutData);
      Real low = OutData(index);
      // If we've specified use_T_floor and we are below Tmin just use Tmin and return
      if (use_T_floor && var < low) {
        return;
      }
      QWData(rho, BrakT[1], OutData);
      Real high =  OutData(index);
      std::stringstream msg;
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("0 ERR (%s): %.4e !<= %.4e !<= %.4e\n", varnames[index], low, var,
             high);
      printf("at rho = %.4e, T_bounds = %.4e, %.4e,\n", rho, BrakT[0], BrakT[1]);
      printf("T_high/low = %.4e, %.4e,\n", T_floor, T_high);
      printf("Init error: %.4e, %.4e\n", a, b);
      msg << "### FATAL ERROR in EquationOfState inversion (TempInvert)"
          << std::endl << "Root not bracketed" << std::endl;
      ATHENA_ERROR(msg);
    }
    // if we are outside brackets use bisection method for a step
    if ((GuessTemp <= BrakT[0]) || (GuessTemp >= BrakT[1])) {
      //GuessTemp = 0.5 * (BrakT[0] + BrakT[1]);
      GuessTemp = std::sqrt(BrakT[0] * BrakT[1]);
#ifdef MYDEBUG1
      mode = 1;
#endif
    }
    QWData(rho, GuessTemp, OutData);
    error = OutData(index) * InvVar - 1.0;
    //error = 1.0 - OutData(index) * InvVar;
#ifdef MYDEBUG1
    printf("%04d [%.4g, %.4g, %.4g]; %.4g| %d\n", 1000 - nlim, BrakT[0], GuessTemp,
           BrakT[1], error, mode);
#endif
    // update bracketing values
    if (error < 0) {
      BrakT[0] = GuessTemp;
      BrakVal[0] = error;
    } else {
      BrakT[1] = GuessTemp;
      BrakVal[1] = error;
    }
    if (BrakT[1] <= BrakT[0]) {
      QWData(rho, BrakT[0], OutData);
      Real low = OutData(index);
      // If we've specified use_T_floor and we are below Tmin just use Tmin and return
      if (use_T_floor && var < low) {
        return;
      }
      QWData(rho, BrakT[1], OutData);
      Real high = OutData(index);
      std::stringstream msg;
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("1 ERR (%s): %.4e !<= %.4e !<= %.4e\n", varnames[index], low, var,
             high);
      printf("at rho = %.4e, T_bounds = %.4e, %.4e,\n", rho, BrakT[0], BrakT[1]);
      printf("T_high/low = %.4e, %.4e,\n", T_floor, T_high);
      printf("Init error: %.4e, %.4e\n", a, b);
      msg << "### FATAL ERROR in EquationOfState inversion (TempInvert)"
          << std::endl << "Root not bracketed" << std::endl;
      ATHENA_ERROR(msg);
    }
    if (std::abs(error) > 1e-2) {
      // secant method step
      delta = error * (GuessTemp - LastTemp) / (error - LastErr);
      LastTemp = GuessTemp;
      GuessTemp -= delta;
      LastErr = error;
#ifdef MYDEBUG1
      mode = 2;
#endif
    } else {
      // Newton–Raphson step
      delta = var * error / OutData(index + 1);
      GuessTemp -= delta;
#ifdef MYDEBUG1
      mode = 3;
#endif
    }
    if (nlim-- < 0) {
      QWData(rho, BrakT[0], OutData);
      Real low = OutData(index);
      QWData(rho, BrakT[1], OutData);
      Real high = OutData(index);
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("2 ERR (%s): |%.4e|, |%.4e| > %.4e; %d iterations\n", varnames[index],
             low * InvVar - 1.0, high * InvVar - 1.0, prec, nmax);
             //1.0 - low * InvVar, 1.0 - high * InvVar, prec, nmax);
      printf("%.4e <= %.4e <= %.4e\n", low, var, high);
      printf("at rho = %.4e, T_bounds = %.4e, %.4e,\n", rho, BrakT[0], BrakT[1]);
      printf("T_high/low = %.4e, %.4e,\n", T_floor, T_high);
      printf("Init error: %.4e, %.4e\n", a, b);
      std::stringstream msg;
      msg << "### FATAL ERROR in EquationOfState inversion (TempInvert)"
          << std::endl << "Cannot converge" << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  if (OutData(5) < T_floor) {
    std::stringstream msg;
    msg << "### FATAL ERROR in EquationOfState inversion (TempInvert)"
        << std::endl << "Cannot converge. Recovered T off table." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//  \brief Return gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {
  TempInvert(rho, LastTemp, egas, EOS::iE, EosData);
  LastTemp = EosData(EOS::iT);
  return EosData(EOS::iP);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//  \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {
  TempInvert(rho, LastTemp, pres, EOS::iP, EosData);
  LastTemp = EosData(EOS::iT);
  return EosData(EOS::iE);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//  \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {
  TempInvert(rho, LastTemp, pres, EOS::iP, EosData);
  LastTemp = EosData(EOS::iT);
  return EosData(EOS::iAsq);
}

void EquationOfState::SevenFromRhoT(Real rho, Real T, AthenaArray<Real> &out) {
  QWData(rho, T, out);
}

Real EquationOfState::TFromRhoP(Real rho, Real pres) {
  TempInvert(rho, LastTemp, pres, EOS::iP, EosData);
  LastTemp = EosData(EOS::iT);
  return LastTemp;
}

Real EquationOfState::TFromRhoEgas(Real rho, Real egas) {
  TempInvert(rho, LastTemp, egas, EOS::iE, EosData);
  LastTemp = EosData(EOS::iT);
  return LastTemp;
}

Real EquationOfState::PresFromRhoT(Real rho, Real T) {
  QWData(rho, T, EosData);
  return EosData(EOS::iP);
}

Real EquationOfState::GetEgasFloor(Real rho) {
  QWData(rho, T_floor, EosData);
  return EosData(EOS::iE);
}

Real EquationOfState::GetPresFloor(Real rho) {
  QWData(rho, T_floor, EosData);
  return EosData(EOS::iP);
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::InitEosConstants(ParameterInput* pin)
//! \brief Initialize constants for EOS
void EquationOfState::InitEosConstants(ParameterInput* pin) {
  EosData.NewAthenaArray(EOS::N);
  prec = pin->GetOrAddReal("hydro", "InversionPrecision", prec);
  Ye = pin->GetReal("problem","Ye");
  T_floor = pin->GetOrAddReal("hydro", "T_floor", float_eps);
  T_ceil = pin->GetOrAddReal("hydro", "T_ceil", float_max);
  use_T_floor = pin->GetOrAddBoolean("hydro", "use_T_floor", true);
  LastTemp = pin->GetOrAddReal("hydro", "T0", std::sqrt(T_floor * T_ceil));
  nmax = pin->GetOrAddInteger("hydro", "nmax", 10000);
  // if root processor and zeroth block
  if ((Globals::my_rank == 0) && (pmy_block_->lid == 0)){
    std::cout<<Ye<<" Ye\n";
  }
  #ifdef MYDEBUG
    Real rho, temp;
    std::cout << "Input fluid parameters and retrieve EOS parameters." << '\n'
              << "Non-positive inputs will exit loop." << '\n';
    std::cout << "Input density (mass/volume): ";
    std::cin >> rho;
    std::cout << "Input temperature (K): ";
    std::cin >> temp;

    while (rho > 0 && std::isfinite(rho) && temp > 0 && std::isfinite(temp)) {
      printf("d, t: %.16e, %.16e\n", rho, temp);
      QWData(rho, temp, EosData);
      Real p = EosData(EOS::iP);
      Real e = EosData(EOS::iE);
      Real a2 = EosData(EOS::iAsq);
      Real P_tejas = compare::P_of_rho_T(rho, temp);
      Real e_tejas = compare::e_of_rho_T(rho, temp);
      Real Asq_tejas = compare::asq(rho, temp, P_tejas);
      printf("e(d, T)          , p(d, T)         , ASq(d, T)\n");
      printf("%.16e, %.16e, %.16e\n", e, p, a2);
      printf("e_matt/e_tejas-1, P_matt/P_tejas-1, Asq_matt/Asq_tejas-1\n");
      printf("%.16e, %.16e, %.16e\n", e/e_tejas-1.0, p/P_tejas-1.0, a2/Asq_tejas-1.0);
      Real Te = TFromRhoEgas(rho, e);
      Real Tp = TFromRhoP(rho, p);
      printf("T(d, e)/T-1, T(d, p)/T-1\n");
      printf("%.16e, %.16e\n", Te/temp - 1.0, Tp/temp - 1.0);
      std::cout << "Input density (mass/volume): ";
      std::cin >> rho;
      std::cout << "Input temperature (K): ";
      std::cin >> temp;
    }
    std::cout << std::endl;
  #endif
  return;
}
