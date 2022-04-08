//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file ideal.cpp
//! \brief implements ideal EOS in general EOS framework, mostly for debuging
//======================================================================================

// C headers

// C++ headers
#include <sstream>
#include <string>

// Athena++ headers
#include "../eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//! \brief Return gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {
  return (gamma_ - 1.) * egas;
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//! \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {
  return pres / (gamma_ - 1.);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//! \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {
  return gamma_ * pres / rho;
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::InitEosConstants(ParameterInput* pin)
//! \brief Initialize constants for EOS
void EquationOfState::InitEosConstants(ParameterInput *pin) {
  return;
}

void EquationOfState::SevenFromRhoT(Real rho, Real T, AthenaArray<Real> &out) {
  std::stringstream msg;
  msg << "### FATAL ERROR in EquationOfState::SevenFromRhoT" << std::endl
      << "Function should not be called with current configuration." << std::endl;
  ATHENA_ERROR(msg);
}
Real EquationOfState::TFromRhoP(Real rho, Real pres) {
  std::stringstream msg;
  msg << "### FATAL ERROR in EquationOfState::TFromRhoP" << std::endl
      << "Function should not be called with current configuration." << std::endl;
  ATHENA_ERROR(msg);
  return -1.0;
}

Real EquationOfState::TFromRhoEgas(Real rho, Real egas) {
  std::stringstream msg;
  msg << "### FATAL ERROR in EquationOfState::TFromRhoEgas" << std::endl
      << "Function should not be called with current configuration." << std::endl;
  ATHENA_ERROR(msg);
  return -1.0;
}
