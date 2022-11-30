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
static Real r_0, inv_r2, rho_0, rho_f, v_f, p_f, mu, Ye, Na, p_0, p0, Tc, Th, tau, T_0, dpdd_0, B_0, alpha, T_floor, Mdot, Factor;
static Real Coeff_nu_0, Coeff_nubar_0, t_L_0, t_L_1, t_coeff, dCoeff_nu, dCoeff_nubar, eps_nue, eps_nueb, L_nue, L_nueb;
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
// zeroQ_temp parameters:
static Real temp_g, dT0, tol, mod, maxc, eps;

static const double kbol      = 1.3806504*std::pow(10,-16);   // erg / K
static const double c         = 2.99792458*std::pow(10,10);   // cm / s
static const double hbar      = 1.05457266*std::pow(10,-27);  // erg s
static const double mn        = 1.6749286*std::pow(10,-24);   // g
static const double electron  = 1.6021773*std::pow(10,-6);    // convert MeV to erg
static const double melectron = 9.1093897*std::pow(10,-28);   // g
static const double delta     = 1.2935;              // MeV
static const double pi        = 3.1415926535897932384626433832795;

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

Real fermi_approx(Real n, Real eta){
  if (n==0) {
    Real fermi = log(1.0 + exp(eta));
    return fermi;
  } else if (n==1) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,2) / 2.0 + 1.6449341;
    Real ff = a - std::pow(a,2) / 4.0 + std::pow(a,3) / 9.0 - std::pow(a,4) / 16.0 + std::pow(a,5) / 25.0 - std::pow(a,6) / 36.0 + std::pow(a,7) / 49.0;
    if (eta < 0) {
        Real fermi = ff;
        return fermi;
    } else if (eta == 0) {
        Real fermi = s - ff;
        return fermi;
    } else {
        Real fermi = s - ff;
        return fermi;
    }
  } else if (n==2) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,3) / 3.0 + 3.2898681 * eta;
    Real ff = 2.0 * (a - std::pow(a,2) / 8.0 + std::pow(a,3) / 27.0 - std::pow(a,4) / 64.0 + std::pow(a,5) / 125.0 - std::pow(a,6) / 216.0);
    if (eta<0) {
        Real fermi = ff;
        return fermi;
    } else if (eta==0) {
        Real fermi = s + ff;
        return fermi;
    } else {
        Real fermi = s + ff;
        return fermi;
    }
  } else if (n==3) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,4) / 4.0 + 4.9348022 * std::pow(eta,2) + 11.351273;
    Real ff = 6.0 * (a - std::pow(a,2) / 16.0 + std::pow(a,3) / 81.0 - std::pow(a,4) / 256.0);
    if (eta<0) {
        Real fermi = ff;
        return fermi;
    } else if (eta==0) {
        Real fermi = s - ff;
        return fermi;
    } else {
        Real fermi = s - ff;
        return fermi;
    }
  } else if (n==4) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,5) / 5.0 + 6.5797363 * std::pow(eta,3) + 45.457576 * eta;
    Real ff = 24.0 * (a - std::pow(a,2) / 32.0 + std::pow(a,3) / 243.0);
    if (eta<0) {
        Real fermi = ff;
        return fermi;
    } else if (eta==0) {
        Real fermi = s + ff;
        return fermi;
    } else {
        Real fermi = s + ff;
        return fermi;
    }
  } else if (n==5) {
    Real a  = exp(-1.0 * fabs(eta));
    Real s  = std::pow(eta,6) / 6.0 + 8.2246703 * std::pow(eta,4) + 113.64394 * std::pow(eta,2) + 236.53226;
    Real ff = 120.0 * (a - std::pow(a,2) / 64.0 + std::pow(a,3) / 729.0);
    if(eta<0) {
        Real fermi = ff;
        return fermi;
    } else if(eta==0) {
        Real fermi = s - ff;
        return fermi;
    } else {
        Real fermi = s - ff;
        return fermi;
    }
  } else {
      std::cout << " (accretion_scheckqdot_QWeos.cpp) \n n for fermi_approx wasnt between 0 and 5. \n Something might not be correct";
      return 0.0;
  }
}

static const double fermi2      = fermi_approx(2,0.0);
static const double fermi4      = fermi_approx(4,0.0);
static const double fermi5      = fermi_approx(5,0.0);
static const double fermi3      = fermi_approx(3,0.0);
static const double fermifactor = fermi2 * fermi5 / ( fermi3 * fermi4 );

Real QWEta(Real rho, Real T, Real Ye) {
  // Returns eta = mu_e / kbT

  // constants
  Real third = 1.0 / 3.0;
  Real c = 2.99792458e10;
  Real k = 1.380649e-16;
  Real mn = 1.6726e-24;
  Real hbar = 6.62607015e-27/(2.0*PI);
  Real c3 = std::pow(k/(hbar*c),3);
  Real eta_den_const = std::pow(6.0,2.0*third);
  Real root3 = std::sqrt(3.0);
  Real eta_den_A = std::pow(2.0,third) / eta_den_const;
  Real eta_den_B = 2.0*std::pow(3.0,third)/eta_den_const;

  Real vol = mn/rho;
  Real T3 = T*T*T;
  Real T4 = T*T3;
  Real a = c3*std::pow(T,3)*vol*(PI/3.0);
  Real a2 = SQR(a);
  Real a4 = SQR(a2);
  Real a6 = a2 * a4;
  Real y2 = SQR(Ye);
  Real b = std::sqrt(4.0*a6+27.0*a4*y2);
  Real term = std::pow(9.0*a2*Ye+root3*b, third);
  Real eta_by_pi = (eta_den_A)*term/a - (eta_den_B)*a/term; // actually eta/pi
  Real eta = eta_by_pi * PI;
  return eta;
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

// Ye sink/source function
Real YeSourceSink(Real temp, Real ye, Real x, Real Eta) {
  // Returns v*dYe/dr (equation 7 in Qian & Woosley 1996) in s^-1 units
  // temp is in MeV
  // Define constants
  Real alpha    = 1.26; // Coupling Coefficient
  Real G_F      = std::pow(298.2*std::pow(10,3),-2.0);  // Fermi Constant in MeV^-2
  Real e_nue    = 1.2 * eps_nue; // Electron neutrino energy defined in QW1996 appendix in MeV
  Real e_nueb   = 1.2 * eps_nueb; // Electron antineutrino energy defined in QW1996 appendix in MeV
  Real Delta    = 1.293; // neutron-proton mass difference in MeV
  Real hbar     = 6.582119569e-22; // Planck's constant in MeV s
  Real c        = 2.99792458e10; // Light speed in cm/s
  Real erg2MeV  = 6.24151e5; // convert erg to MeV
  Real lum      = std::pow(10,51);


  // Forward reaction rates
  Real lambda_nue_n  = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*L_nue*erg2MeV*lum/(r_0*r_0)
                      *(e_nue+2.0*Delta+1.2*Delta*Delta/e_nue)*(1.0-x); // s^-1 units
  Real lambda_nueb_p = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*L_nueb*erg2MeV*lum/(r_0*r_0)
                      *(e_nueb-2.0*Delta+1.2*Delta*Delta/e_nueb)*(1.0-x); // s^-1 units

  // Reverse reaction rates
  Real f4_0 = fermi_approx(4.0,0.0);
  Real f4_eta = fermi_approx(4.0,Eta);
  Real f4_negeta = fermi_approx(4.0,-1.0*Eta);
  Real lambda_ele_p = 0.448*std::pow(temp,5.0)*f4_eta/f4_0; // s^-1 units (e- p)
  Real lambda_pos_n = lambda_ele_p*f4_negeta/f4_0; // s^-1 units (e+ n)

//  std::cout << "T             = " << temp << " MeV" << std::endl;
//  std::cout << "Ye            = " << ye << std::endl;
//  std::cout << "f4_0          = " << f4_0 << std::endl;
//  std::cout << "f4_eta        = " << f4_eta << std::endl;
//  std::cout << "f4_negeta     = " << f4_negeta << std::endl;
//  std::cout << "lambda_nue_n  = " << lambda_nue_n << " s^-1" << std::endl;
//  std::cout << "lambda_nueb_p = " << lambda_nueb_p << " s^-1" << std::endl;
//  std::cout << "lambda_ele_p  = " << lambda_ele_p << " s^-1" << std::endl;
//  std::cout << "lambda_pos_n  = " << lambda_pos_n << " s^-1" << std::endl;
//  std::cout << " " << std::endl;

  // Source - Sink
  Real out = lambda_nue_n+lambda_pos_n-(lambda_nue_n+lambda_pos_n+lambda_nueb_p+lambda_ele_p)*ye;
  return out;
}

Real Ye_f(Real temp, Real ye, Real x){
  // Returns v*dYe/dr (equation 7 in Qian & Woosley 1996) in s^-1 units
  // temp is in MeV
  // Define constants
  Real alpha  = 1.26; // Coupling Coefficient
  Real G_F    = std::pow(298.2*std::pow(10,3),-2.0);  // Fermi Constant in MeV^-2
  Real e_nue  = 1.2 * eps_nue; // Electron neutrino energy defined in QW1996 appendix in MeV
  Real e_nueb = 1.2 * eps_nueb; // Electron antineutrino energy defined in QW1996 appendix in MeV
  Real Delta  = 1.293; // neutron-proton mass difference in MeV
  Real hbar   = 6.582119569e-22; // Planck's constant in MeV s
  Real c      = 2.99792458e10; // Light speed in cm/s
  Real erg2MeV = 6.24151e5; // convert erg to MeV
  Real lum     = std::pow(10,51);

  // Forward reaction rates
  Real lambda_nue_n  = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*L_nue*erg2MeV*lum/(r_0*r_0)
                      *(e_nue+2.0*Delta+1.2*Delta*Delta/e_nue)*(1.0-x); // s^-1 units
  Real lambda_nueb_p = std::pow(hbar*c,2.0)*(1.0+3.0*alpha*alpha)/(2*PI*PI)*G_F*G_F*L_nueb*erg2MeV*lum/(r_0*r_0)
                      *(e_nueb-2.0*Delta+1.2*Delta*Delta/e_nueb)*(1.0-x); // s^-1 units

  // Source - Sink (look at lnue and lnueb, they should be smaller)
  Real out = (lambda_nue_n) / (lambda_nue_n + lambda_nueb_p);
//  std::cout << "alpha = " << alpha << std::endl;
//  std::cout << "G_F   = " << G_F << " MeV^-2" << std::endl;
//  std::cout << "e_nue = 1.2 * eps_nue = " << e_nue << " MeV" << std::endl;
//  std::cout << "e_nueb = 1.2 * eps_nueb = " << e_nueb << " MeV" << std::endl;
//  std::cout << "Delta = " << Delta << " MeV" << std::endl;
//  std::cout << "hbar = " << hbar << " MeV s" << std::endl;
//  std::cout << "c = " << c << " cm/s" << std::endl;
//  std::cout << "Lnue = " << L_nue*erg2MeV*1.0e51 << " MeV/s" << std::endl;
//  std::cout << "Lnueb = " << L_nueb*erg2MeV*1.0e51 << " MeV/s" << std::endl;
//  std::cout << "r_0 = " << r_0 << " cm" << std::endl;
//  std::cout << "lambda_nue_n = " << lambda_nue_n << " s^-1" << std::endl;
//  std::cout << "lambda_nueb_p = " << lambda_nueb_p << " s^-1" << std::endl;
//  std::cout << "Y_e,f = " << out << std::endl;

  return out;
}

Real qdotScheck_H(Real ye, Real r) {
  // returns heating in units of erg s^{-1} g^{-1} units
  const double alpha     = 1.254;
  const double sigmaweak = 1.76*std::pow(10,-44)*(1.0+3.0*std::pow(alpha,2))/(4.0*std::pow(melectron*c*c,2));   // weak cross section cm^2 / erg^2
  Real Lnu               = L_nue;

  Real fnu         = 0.5 * (1.0 + std::pow(1.0 - std::pow(r_0 / r,2),0.5));
  Real sinkconst   = sigmaweak * 0.5 * (Lnu * std::pow(10,51)) / (4.0 * PI * std::pow(r,2) * fnu); // erg^-1 s^-1

  // --- Absorption Term ---
  Real q_abs  = sinkconst / mn * ((1.0 - ye) * (std::pow(eps_nue * electron,2) * fermifactor +
                                   2.0 * delta * electron * eps_nue * electron * std::pow(fermi_approx(2,0) * fermi_approx(4,0),0.5) / fermi_approx(3,0) +
                                   std::pow(delta * electron,2)) + ye * (std::pow(eps_nueb * electron,2) * fermifactor +
                                   3.0 * delta * electron * eps_nueb * electron * std::pow(fermi_approx(2,0) * fermi_approx(4,0),0.5) / fermi_approx(3,0) +
                                   3.0 * std::pow(delta * electron,2))); // erg s^-1 g^-1

//  std::cout << " ";
//  std::cout << "qdotScheck_H\n";
//  std::cout << "sigmaweak = " << sigmaweak << std::endl;
//  std::cout << "fnu       = " << fnu       << std::endl;
//  std::cout << "r_0       = " << r_0       << " cm" << std::endl;
//  std::cout << "r         = " << r         << " cm" << std::endl;
//  std::cout << "sinkconst = " << sinkconst << std::endl;
//  std::cout << "delta       = " << delta << std::endl;
//  std::cout << "fermifactor = " << fermifactor << std::endl;
//  std::cout << "eps_nue     = " << eps_nue << std::endl;
//  std::cout << "eps_nueb    = " << eps_nueb << std::endl;
//  std::cout << "electron    = " << electron << std::endl;
//  std::cout << "q_abs     = " << q_abs << std::endl;
//  std::cout << "Lnu       = " << Lnu << std::endl;
//  std::cout << "Ye (qdotS) = " << ye << std::endl;
////  std::cout << "r  (qdotS) = " << r << std::endl;
//  std::cout << " " << std::endl;


  return q_abs;
}

Real qdotScheck_C(Real temp, Real ye, Real rho, Real r) {
  // temp is in K
  // returns cooling in units of erg s^{-1} g^{-1} units
  const double alpha       = 1.254;
  const double sigmaweak   = 1.76*std::pow(10,-44)*(1.0+3.0*std::pow(alpha,2))/(4.0*std::pow(melectron*c*c,2));   // weak cross section cm^2 / erg^2

  Real fnu         = 0.5 * (1.0 + std::pow(1.0 - std::pow(r_0 / r,2),0.5));
  Real Tprime      = temp * kbol / (hbar * c); // cm^-1
  Real Terg        = temp * kbol; // erg
  Real tfact       = std::pow(27.0 * ye * rho + std::pow(3.0,0.5) * std::pow(4.0 * std::pow(pi,2) * std::pow(mn,2) * std::pow(Tprime,6) + 243.0 * std::pow(ye,2) * std::pow(rho,2),0.5),1.0/3.0); // g^(1/3) cm^-1
  Real eta_e       = std::pow(12.0 * pi,1.0/3.0) * (std::pow(pi * mn,1.0/3.0) * std::pow(tfact,2) - pi * std::pow(Tprime,2) * mn * std::pow(12.0,1.0/3.0)) / (6.0 * std::pow(mn,2.0/3.0) * tfact * Tprime); // unitless
  Real sourceconst = sigmaweak * 0.5 * c * std::pow(Tprime,3) / (std::pow(pi,2) * mn); // erg^-2 s^-1 g^-1

  // --- Emission Term ----
  Real q_em   = sourceconst * (ye * (std::pow(Terg,3) * fermi_approx(5,eta_e - delta * electron / Terg) +
                              2.0 * delta * electron * std::pow(Terg,2) * fermi_approx(4,eta_e - delta * electron / Terg) +
                              std::pow(delta * electron,2) * Terg * fermi_approx(3,eta_e - delta * electron / Terg)) +
                              (1.0 - ye) * (std::pow(Terg,3) * fermi_approx(5,-1.0*eta_e) + 3.0 * delta * electron * std::pow(Terg,2) * fermi_approx(4,-1.0*eta_e) + 
                              3.0 * std::pow(delta * electron,2) * Terg * fermi_approx(3,-1.0*eta_e) + std::pow(delta * electron,3) * fermi_approx(2,-1.0*eta_e))); // erg s^-1 g^-1

//  std::cout << " ";
//  std::cout << "qdotScheck_C\n";
//  std::cout << "sigmaweak   = " << sigmaweak << "\n";
//  std::cout << "fnu         = " << fnu       << "\n";
//  std::cout << "r_0       = " << r_0       << " cm" << std::endl;
//  std::cout << "r         = " << r         << " cm" << std::endl;
//  std::cout << "Tprime      = " << Tprime << "\n";
//  std::cout << "Terg        = " << Terg << "\n";
//  std::cout << "tfact       = " << tfact << "\n";
//  std::cout << "eta_e       = " << eta_e << "\n";
//  std::cout << "sourceconst = " << sourceconst << "\n";
//  std::cout << "q_em        = " << q_em << " erg/g/s\n";
//  std::cout << " ";

  return q_em;
}


// Temperature at which heating=cooling
//inline Real zeroQ_temp(Real x, Real ye, Real time=0.0) {
//  // smoothly transition from 0 at t<=t_L_0 to 1 at t>=t_L_1
//  Real f = (time <= t_L_0) ? 0.0 :
//           ((time >= t_L_1) ? 1.0 : SQR(std::sin((time - t_L_0) * t_coeff)));
//  Real Coeff_nu = Coeff_nu_0 + f * dCoeff_nu;
//  Real Coeff_nubar = Coeff_nubar_0 + f * dCoeff_nubar;
//  return std::pow(1e12*9.65/2.27*((1.-ye)*Coeff_nu+ye*Coeff_nubar)*(1.-x)*inv_r2,
//                  1./6.) / 8.6173e-11;
//}

inline Real zeroQ_temp(Real ye, Real rho, Real r) {
  // Perform NR on qdot to determine T at r=R_PNS such that H=C
  Real T0 = temp_g;
  for(int i=0; i<=maxc; i++) {
    Real DeltaT0 = dT0*T0;
    Real y       = (qdotScheck_H(ye,r) - qdotScheck_C(T0,ye,rho,r));
    Real yprime  = ((qdotScheck_H(ye,r) - qdotScheck_C(T0+DeltaT0,ye,rho,r)) -
                   (qdotScheck_H(ye,r) - qdotScheck_C(T0,ye,rho,r))) / ((T0+DeltaT0) - (T0));
//    std::cout << "qdotScheck_H = " << qdotScheck_H(ye,r) << std::endl;
    if (fabs(yprime) < eps) {
      std::cout << "(zeroQ_temp) yprime is too small!\n";
      return T0;
      break;
    }
    Real T1 = T0 - mod * y / yprime;
    if (fabs(T1-T0) / T1 < tol) {
//      std::cout << "Root found!\n";
//      std::cout << "T1 (zeroQ_T) = " << T1 << " K" << std::endl;
      return T1;
      break;
    }
    T0 = T1;
    if (i==maxc) {
      std::cout << "(zeroQ_temp) maximum iterations reached.\n";
      return T0;
      break;
    }
  }
  return T0;
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
  Mdot = pin->GetReal("problem","mdot");
  Factor = pin->GetOrAddReal("problem","rho0factor",1.0);
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
  ye_const = pin->GetReal("hydro","Ye_f");
  ye_0     = pin->GetReal("hydro","Ye_0");
  // zeroQ parameters
  temp_g = pin->GetReal("problem","Temp_guess");
  dT0    = pin->GetReal("problem","dTemp");
  tol    = pin->GetReal("problem","tolerance");
  mod    = pin->GetReal("problem","modifier");
  maxc   = pin->GetReal("problem","maxcount");
  eps    = pin->GetReal("problem","epsilon");
  // final limonosity/energies
  Real L_nu = pin->GetReal("problem","L_nu");
  Real L_nubar = pin->GetReal("problem","L_nubar");
  L_nue = pin->GetReal("problem","L_nu");
  L_nueb = pin->GetReal("problem","L_nubar");
  Real eps_nu = pin->GetReal("problem","eps_nu");
  Real eps_nubar = pin->GetReal("problem","eps_nubar");
  eps_nue = pin->GetReal("problem","eps_nu");
  eps_nueb = pin->GetReal("problem","eps_nubar");
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
    int col_Ye = pin->GetInteger("problem", "col_Ye");

    // Prepare arrays to hold profile
    AllocateRealUserMeshDataField(7);
    ruser_mesh_data[0].NewAthenaArray(rows);
    ruser_mesh_data[1].NewAthenaArray(rows);
    ruser_mesh_data[2].NewAthenaArray(rows);
    ruser_mesh_data[3].NewAthenaArray(rows);
    ruser_mesh_data[4].NewAthenaArray(rows);
    AthenaArray<Real>& r_in{ruser_mesh_data[0]};
    AthenaArray<Real>& rho_in{ruser_mesh_data[1]};
    AthenaArray<Real>& v_in{ruser_mesh_data[2]};
    AthenaArray<Real>& T_in{ruser_mesh_data[3]};
    AthenaArray<Real>& Ye_in{ruser_mesh_data[4]};

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
      Ye_in(n) = s_vals[col_Ye+1];
      if (Globals::my_rank == 0) {
       // std::cout<<r_in(n)<<" ,"<<rho_in(n)<<" , "<<v_in(n)<<" , "<<T_in(n);
       // std::cout<<"\n";
      }
    }

  }
  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  Real ne_init = ye_0*rho_0; //pscalars->s(ye_index,0,0,0) = ye_0*rho_0;
  Real ye_init = ye_0; //pscalars->r(ye_index,0,0,0) = pscalars->s(ye_index,0,0,0)/rho_0;
  edens0 = &ne_init; //&pscalars->s(ye_index,0,0,0);
  Ye0    = &ye_init; //&pscalars->r(ye_index,0,0,0);
  Real r = pcoord->x1v(0);
//  T_0    = zeroQ_temp(0.0, ye_0);
  T_0    = zeroQ_temp(ye_0, rho_0, r_0);
  p_0    = peos->PresFromRhoT(rho_0, T_0, edens0);
  std::cout << "T_0 = " << T_0 << " K" << std::endl;
  std::cout << "p_0 = " << p_0 << " erg/cm^3" << std::endl;
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
    AthenaArray<Real>& Ye_in{pmy_mesh->ruser_mesh_data[4]};


    for (int k=ks; k<=ke; k++) {
      // Real phi = pcoord->x3v(k);
      for (int j=js; j<=je; j++) {
        Real theta = pcoord->x2v(j);
     //	std::cout<< "x2v(j) = " << pcoord->x2v(j) << "\n";
        for (int i=is; i<=ie; i++) {
          Real r = pcoord->x1v(i);
          Real r0 = 5e6;
          Real rho, v, temp, ye;

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
          ye=Ye_in(index)+(r-r_in(index))*(Ye_in(index+1)-Ye_in(index))
               /(r_in(index+1)-r_in(index));

          phydro->u(IDN,k,j,i) = rho;
          phydro->u(IM1,k,j,i) = v * rho;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          pscalars->s(ye_index,k,j,i) = ye*rho;
          pscalars->r(ye_index,k,j,i) = ye;
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
//  Real sumye = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real erg2MeV  = 6.24151e5;
        Real r = pmb->pcoord->x1v(i);
        Real p = prim(IPR,k,j,i);
        Real rho = prim(IDN,k,j,i);

        edens = &cons_scalar(ye_index,k,j,i);
        Real temp = pmb->peos->TFromRhoP(rho, p, edens) * 8.6173e-11; //T in MeV
        Real x = std::sqrt(1.0-(r_0*r_0)/(r*r));
        Real Ye_i = cons_scalar(ye_index,k,j,i)/prim(IDN,k,j,i);
        Real eta = QWEta(rho,temp/8.6173e-11,Ye_i);
        //std::cout << "eta = " << eta << std::endl;
        Real vdYedr = YeSourceSink(temp,Ye_i,x,eta);
        Real dYe = dt*vdYedr;
//        sumye += dYe;
//        std::cout << "sumye (step) = " << sumye << std::endl;
//        std::cout << "dt               = " << dt << std::endl;
//        std::cout << "dYe = dt*vdYe/dr = " << dYe << std::endl;
//        std::cout << " " << std::endl;
        cons_scalar(ye_index,k,j,i) += dYe*prim(IDN,k,j,i);
//        if (i==pmb->is) {
//          std::cout << "Ye                = " << Ye_i << std::endl;
//          std::cout << "Ye (scalar field) = " << prim_scalar(ye_index,k,j,i) << std::endl;
//        }

        // Real testingthing = Ye_f(temp, Ye_i, x);
        //Real elecfrac = cons_scalar(ye_index,k,j,i)/prim(IDN,k,j,i);

        //Real qdot = qdotQW(temp, x, elecfrac, time); //MeV s^{-1} g^{-1} units
        //Real qdot = qdotQW(temp, x, Ye_i, time)/erg2MeV;  //erg s^-1 g^-1 units
        Real qdot = qdotScheck_H(Ye_i,r) - qdotScheck_C(temp/8.6173e-11,Ye_i,rho,r); // erg/s/g units
        //std::cout << "qdot = " << qdot << " erg/s/g" << std::endl;
        Real de = dt*prim(IDN,k,j,i) * qdot; // * 1.6022e-6; //removed tanh
        cons(IEN,k,j,i) += de;
      }
    }
  }
//  std::cout << " " << std::endl;
//  std::cout << "Sum Ye = " << sumye << std::endl;
//  std::cout << " " << std::endl;
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
        //Real r = pco->x1v(is);
	Real r = pco->x1v(is-i);
//	std::cout << "r_0           = " << r_0 << " cm" << std::endl;
//	std::cout << "r             = " << r << " cm " << std::endl;
//	std::cout << "exponent term = " << std::exp((r_0 - r) * mu / (dpdd_0 * r * r_0)) << std::endl;
//        prim(IDN,k,j,is-i) = rho_0;
        prim(IDN,k,j,is-i) = Factor * rho_0 * std::exp((r_0 - r) * mu / (dpdd_0 * r * r_0));
//	std::cout << "rho_ghost = " << prim(IDN,k,j,is-i) << " g/cm^3" << std::endl;
//      prim(IVX,k,j,is-i) = std::max(prim(IVX,k,j,is), 0.0);
        prim(IVX,k,j,is-i) = (-1.0 * Mdot * (2.0e33)) / (4.0 * PI * std::pow(r_0,2) * rho_0);
//        prim(IVX,k,j,is-i) = (-1.0 * Mdot * (2.0e33)) / (4.0 * PI * std::pow(pco->x1v(is),2) * prim(IDN,k,j,is));
//        prim(IVX,k,j,is-i) = prim(IVX,k,j,is);
        prim(IVY,k,j,is-i) = 0.0;
        prim(IVZ,k,j,is-i) = 0.0;
        // Setting inner ghost zone Ye equal to first active zone value
        pmb->pscalars->r(ye_index,k,j,is-i) = pmb->pscalars->r(ye_index,k,j,is);
//	Real ye=pmb->pscalars->r(ye_index,k,j,is-i);
//	Real rho=prim(IDN,k,j,is-i);
//        std::cout << "--------------INFLOWINNERX1--------------" << std::endl;
//	Real T_ghost = zeroQ_temp(ye,rho,r);
//	std::cout << "ye      = " << ye << std::endl;
//	std::cout << "rho     = " << rho << " g/cm^3" << std::endl;
//	std::cout << "T_ghost = " << T_ghost << " K" << std::endl;
//        std::cout << "Inner BC" << std::endl;
//        std::cout << "is-i = " << is-i << std::endl;
//        std::cout << "i = " << i << std::endl;
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
//        std::cout << "Outer BC" << std::endl;
//        std::cout << "ie+i = " << ie+i << std::endl;
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
//        Real qdot = qdotQW(temp, x, pscalars->r(ye_index,k,j,i), pmy_mesh->time) * 1.6022e-6; // in ergs/s/g
//        Real qdot = qdotScheck_H(pscalars->r(ye_index,k,j,i),r) -
//                    qdotScheck_C(temp/8.6173e-11,pscalars->r(ye_index,k,j,i),phydro->u(IDN,k,j,i),r);
//        Real yef  = Ye_f(temp, pscalars->r(ye_index,k,j,i), x);
        Real eta = QWEta(phydro->u(IDN<k,j,i),temp/8.6173e-11,pscalars->r(ye_index,k,j,i));
        user_out_var(1,k,j,i) = qdotScheck_H(pscalars->r(ye_index,k,j,i),r) - qdotScheck_C(temp/8.6173e-11,pscalars->r(ye_index,k,j,i),phydro->u(IDN,k,j,i),r);
        user_out_var(2,k,j,i) = std::sqrt(peos->AsqFromRhoP(phydro->w(IDN,k,j,i),
                                                            phydro->w(IPR,k,j,i),efrac));
//        user_out_var(2,k,j,i) = qdotScheck_C(temp/8.6173e-11,pscalars->r(ye_index,k,j,i),phydro->u(IDN,k,j,i),r);
      }
    }
  }
}
