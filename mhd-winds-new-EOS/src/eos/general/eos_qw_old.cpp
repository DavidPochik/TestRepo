//#define MYDEBUG1
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
  int nmax;
  AthenaArray<Real> EosData;

  const Real pi= 3.14159265;
  const Real c = 3e10;
  const Real k = 1.380649e-16;
  const Real mn = 1.6726e-24;
  const Real hbar = 1.0546e-27;
  const Real Tmin= 1e3;
  const Real Tmax= 1e15;
  const Real c3 = std::pow(k/(hbar*c),3);
  const Real con3 = (11.0*pi*pi/180.0)*c3*k;
}


//returns \eta/pi in terms of Ye. See eqn 6 of QW
Real eta_by_pi(Real rho, Real T){
  
  Real a= std::pow(k*T/(hbar*c),3)*(mn/rho)*(pi/3.0);
  Real term= std::pow(9.0*a*a*Ye+1.7321*std::pow((4.0*std::pow(a,6.0)+27.0*std::pow(a,4.0)*Ye*Ye),0.5),1.0/3.0);
  Real r1=0.38157*term/a - 0.87358*a/term;
  return r1;
}

//returns (d \eta/pi)/dT in terms of Ye. See eqn 6 of QW
Real der_etabypi(Real rho, Real T){
 
  Real a= std::pow(k*T/(hbar*c),3)*(mn/rho)*(pi/3.0);
  Real da_dt= 3.0*T*T*std::pow(k/(hbar*c),3)*(mn/rho)*(pi/3.0);
  Real term= std::pow(9.0*a*a*Ye+1.7321*std::pow((4.0*std::pow(a,6.0)+27.0*std::pow(a,4.0)*Ye*Ye),0.5),1.0/3.0);
  Real dterm_dt= (1.0/3.0)*std::pow(9.0*a*a*Ye+1.7321*std::pow((4.0*std::pow(a,6.0)+27.0*std::pow(a,4.0)*Ye*Ye),0.5),-2.0/3.0)*(18.0*a*Ye*da_dt + 1.7321*0.5*std::pow((4.0*std::pow(a,6.0)+27.0*std::pow(a,4.0)*Ye*Ye),-0.5)*(24.0*std::pow(a,5.0)*da_dt+27.0*4.0*a*a*a*Ye*Ye*da_dt));
  return 0.38157*(dterm_dt/a - term*da_dt/(a*a)) - 0.87358*(da_dt/term -a*dterm_dt/(term*term));
  

}
//! \brief compute gas pressure, see eqn 4 of QW
Real P_of_rho_T(Real rho, Real T) {
  Real eta = eta_by_pi(rho,T);
  Real con= (11.0*pi*pi/180.0)*std::pow(k/(hbar*c),3)*k*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0);
  return con*std::pow(T,4.0)+ rho*k*T/mn;
}

//----------------------------------------------------------------------------------------
//! \brief compute internal energy density, see eqn 5 of QW
Real e_of_rho_T(Real rho, Real T) {
  Real eta = eta_by_pi(rho,T);
  Real constant= (11.0*pi*pi/60.0)*std::pow(k/(hbar*c),3)*k*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0);
  Real sp_energy= constant*std::pow(T,4.0)/rho + 1.5*k*T/mn; //energy/g
  return rho*sp_energy;
}



//temp derivative of the function
Real func_derT(Real rho,Real T) {
  Real eta = eta_by_pi(rho,T);

  Real der_eta= der_etabypi(rho,T);
  Real con3= (11.0*pi*pi/180.0)*std::pow(k/(hbar*c),3)*k;
  Real val2= 4.0*con3*std::pow(T,3.0)*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0)+con3*std::pow(T,4.0)*(30.0*2.0*eta*der_eta/11.0 + 15.0*4.0*std::pow(eta,3)*der_eta/11.0)+ rho*k/mn;
  return val2;
}


Real func_derT_egas(Real rho,Real T) {
  Real eta = eta_by_pi(rho,T);
  
  Real der_eta= der_etabypi(rho,T);
  Real con3= (11.0*pi*pi/60.0)*std::pow(k/(hbar*c),3)*k;
  Real val2= (4.0*con3*std::pow(T,3.0)*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0)+con3*std::pow(T,4.0)*(30.0*2.0*eta*der_eta/11.0 + 15.0*4.0*std::pow(eta,3)*der_eta/11.0))+ 1.5*rho*k/mn;
  return val2;
}

//----------------------------------------------------------------------------------------
//! \brief compute adiabatic sound speed squared
Real asq(Real rho, Real T) {
  Real eta = eta_by_pi(rho,T);
  Real ctsq=k*T/mn;
  Real der_eta= der_etabypi(rho,T);
  Real con1= (11.0*pi*pi/60.0)*std::pow(k/(hbar*c),3)*k*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0);
  Real D=(T/rho)*func_derT(rho,T);
  Real con3= (11.0*pi*pi/60.0)*std::pow(k/(hbar*c),3)*k;
  Real cv= (1.0/rho)*(4.0*con3*std::pow(T,3.0)*(1.0+30.0*eta*eta/11.0 + 15.0*std::pow(eta,4)/11.0)+con3*std::pow(T,4.0)*(30.0*2.0*eta*der_eta/11.0 + 15.0*4.0*std::pow(eta,3)*der_eta/11.0))+ 1.5*k/mn;
  return ctsq+D*D/(cv*T);
}


// EOS data indicies
namespace EOS {
  enum EosIndex {iE=0, idEdT=1, iP=2, idPdT=3, iAsq=4, iT=5, N=6};
}



void QWData(Real rho, Real T, AthenaArray<Real> &OutData){

  Real P = P_of_rho_T(rho,T);
  Real e = e_of_rho_T(rho,T);

  Real der_P= func_derT(rho,T);
  Real der_e= func_derT_egas(rho,T);
  Real vsq= asq(rho,T);

  
  OutData(0) = e;
  OutData(1) = der_e;
  OutData(2) = P;
  OutData(3) = der_P;
  OutData(4) = vsq;
  OutData(5) = T;
}



// index = 0 for internal energy; index = 2 for pressure; var = int energy or pressure
void TempInvert(Real rho, Real GuessTemp, Real var, const int index,
                AthenaArray<Real> &OutData) {
  if ((std::isnan(var)) || (std::isnan(rho)) || (std::isnan(GuessTemp))) {
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("ERR (%s): %.4e, rho: %.4e,  Temp: %.4e \n", varnames[index],var,rho, GuessTemp);
      
      std::stringstream msg;
      msg <<"Nan in root find (var) \n";
      ATHENA_ERROR(msg);
    }
  if ((var<=0.0) || (rho<=0.0)) {
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("ERR (%s): %.4e, rho: %.4e,  Temp: %.4e \n", varnames[index],var,rho, GuessTemp);
      std::cout<<eta_by_pi(rho,GuessTemp)<<"    "<<e_of_rho_T(rho,GuessTemp)<<"\n";
      std::stringstream msg;
      msg <<"Negative var\n";
      ATHENA_ERROR(msg);
  }

  Real BrakT[] = {T_floor,T_ceil};
  Real BrakVal[] = {0, 0};
  Real InvVar = 1.0 / var;
  Real error = float_max;
  int nlim = nmax;
  Real LastTemp = BrakT[0];
  QWData(rho, LastTemp, OutData);
  BrakVal[0] = OutData(index) * InvVar - 1.0;
  Real LastErr = BrakVal[0];
  Real delta;
#ifdef MYDEBUG1
  //printf("%d: %.16e, %.16e\n", index, var, rho);
  int mode = 0;
#endif
  while (std::abs(error) > prec) {

  if (BrakVal[0] > 0) {
     	QWData(rho, BrakT[0], OutData);
      	Real low = OutData(index);
      // If we are below Tmin
      	if (var < low) {
	   Real rtemp;
	   if (index == EOS::iP) {
	     rtemp = std::min(std::pow(var / con3, .25), var * mn / (k * rho));
	     QWData(rho, rtemp, OutData);
	     
	   }
	   else if (index == EOS::iE) {
	     rtemp = std::min(std::pow(var / (3.0 * con3), .25), var * mn / (1.5 * k * rho));
	     QWData(rho, rtemp, OutData);
	   }
	   
	   return;
	}
      QWData(rho, BrakT[1], OutData);
      Real high =  OutData(index);
      std::stringstream msg;
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("0 ERR (%s): %.4e !<= %.4e !<= %.4e\n", varnames[index], low, var,
             high);
      printf("at rho = %.4e, T_bounds = %.4e, %.4e,\n", rho, BrakT[0], BrakT[1]);
      
      msg << "### FATAL ERROR in EquationOfState inversion (TempInvert)"
          << std::endl << "Root not bracketed" << std::endl;
      ATHENA_ERROR(msg);

  }
   	  
    // if we are outside brackets use bisection method for a step
   if ((GuessTemp <= BrakT[0]) || (GuessTemp >= BrakT[1])) {
      		if((BrakT[0]<= 0.0) || (BrakT[1]<= 0.0))
	{
	   Real rtemp;
	   if (index == EOS::iP) {
	     rtemp = std::min(std::pow(var / con3, .25), var * mn / (k * rho));
	     QWData(rho, rtemp, OutData);
	     std::cout<<"Var:pres, BrakT<0   "<<rho<<"    "<<var<<"   "<<rtemp<<"\n";
	   }
	   else if (index == EOS::iE) {
	     rtemp = std::min(std::pow(var / (3.0 * con3), .25), var * mn / (1.5 * k * rho));
	     QWData(rho, rtemp, OutData);
	     std::cout<<"Var:egas, BrakT<0   "<<rho<<"    "<<var<<"    "<<rtemp<<"\n";
	   }
	   return;
	}
     
	  
      	GuessTemp = std::sqrt(BrakT[0] * BrakT[1]);
      	if(std::isnan(GuessTemp)){std::cout<<"Nan guessT in brak   "<<BrakT[0]<<"    "<<BrakT[1]<<"\n";}
#ifdef MYDEBUG1
      	mode = 1;
#endif
   }
    	QWData(rho, GuessTemp, OutData);
    	error = OutData(index) * InvVar - 1.0;
#ifdef MYDEBUG1
    	if(GuessTemp<0.0 ||(std::isnan(GuessTemp))){
    		printf("%04d [%.4g, %.4g, %.4g]; %.4g| %d\n", 10000 - nlim, BrakT[0], GuessTemp,
           	BrakT[1], error, mode);
	}
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
      // If we've specified Tfloor and we are below Tmin just use Tmin and return
      		if (var < low) {
		    Real rtemp;
		    if (index == EOS::iP) {
		      rtemp = std::min(std::pow(var / con3, .25), var * mn / (k * rho));
		      QWData(rho, rtemp, OutData);
	    
		    }
		    else if (index == EOS::iE) {
		      rtemp = std::min(std::pow(var / (3.0 * con3), .25), var * mn / (1.5 * k * rho));
		      QWData(rho, rtemp, OutData);
		    }
		    return;
		}
      QWData(rho, BrakT[1], OutData);
      Real high = OutData(index);
      std::stringstream msg;
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      printf("1 ERR (%s): %.4e !<= %.4e !<= %.4e\n", varnames[index], low, var,
             high);
      printf("at rho = %.4e, T_bounds = %.4e, %.4e,\n", rho, BrakT[0], BrakT[1]);
    
      msg << "### FATAL ERROR in EquationOfState inversion (TempInvert)"
          << std::endl << "Root not bracketed" << std::endl;
      ATHENA_ERROR(msg);
	}
  
    if (std::abs(error) > 1e-2) {
      // secant method step
      delta = error * (GuessTemp - LastTemp) / (error - LastErr);
      LastTemp = GuessTemp;
      GuessTemp -= delta;
      
      //if(std::isnan(GuessTemp)){std::cout<<"Nan in secant   "<<rho<<"     "<<delta<<"    "<<error-LastErr<<"    "<<GuessTemp<<"    "<<LastTemp<<"    "<<error<<"\n";}

      LastErr = error;
#ifdef MYDEBUG1
      mode = 2;
#endif
    } else {
      // Newton–Raphson step
      delta = var * error / OutData(index + 1);
      GuessTemp -= delta;
      //if(std::isnan(GuessTemp)){std::cout<<"Nan in NR  "<<OutData(index+1)<<"    "<<delta<<"\n";}

#ifdef MYDEBUG1
      mode = 3;
#endif
    }
    if (nlim-- < 0) {
      	   Real rtemp;
	   if (index == EOS::iP) {
	     rtemp = std::min(std::pow(var / con3, .25), var * mn / (k * rho));
	     QWData(rho, rtemp, OutData);
	     //std::cout<<"Var:pres, nlim reached   "<<rho<<"    "<<var<<"\n";
	   }
	   else if (index == EOS::iE) {
	     rtemp = std::min(std::pow(var / (3.0 * con3), .25), var * mn / (1.5 * k * rho));
	     QWData(rho, rtemp, OutData);
	     //std::cout<<"Var:egas, nlim reached   "<<rho<<"    "<<var<<"\n";
	   }
	   return;
    }
  }

  if (OutData(5) < Tmin) {
           Real rtemp;
	   if (index == EOS::iP) {
	     rtemp = std::min(std::pow(var / con3, .25), var * mn / (k * rho));
	     QWData(rho, rtemp, OutData);
	     //std::cout<<"Var:pres, T<Tmin   "<<rho<<"    "<<var<<"\n";
	   }
	   else if (index == EOS::iE) {
	     rtemp = std::min(std::pow(var / (3.0 * con3), .25), var * mn / (1.5 * k * rho));
	     QWData(rho, rtemp, OutData);
	     //std::cout<<"Var:egas, T<Tmin   "<<rho<<"    "<<var<<"\n";
	   }
     
	  // std::cout<<rho<<"   Return_temp < Tmin encountered, returning Tmin\n";
	   return;
  }
  if (std::isnan(OutData(5))) {
      const char *varnames[] = {"e_int", "de_int/dT", "P_gas", "dP/dT", "a^2", "T"};
      //printf("ERR (%s): %.4e, rho: %.4e \n", varnames[index],var,rho);
      Real rtemp;
      if (index == EOS::iP) {
	     rtemp = std::min(std::pow(var / con3, .25), var * mn / (k * rho));
	     QWData(rho, rtemp, OutData);
	     //std::cout<<"Var:pres, T nan   "<<rho<<"    "<<var<<"    "<<rtemp<<"\n";
	   }
       else if (index == EOS::iE) {
	     rtemp = std::min(std::pow(var / (3.0 * con3), .25), var * mn / (1.5 * k * rho));
	     QWData(rho, rtemp, OutData);
	     //std::cout<<"Var:egas, T nan   "<<rho<<"    "<<var<<"    "<<rtemp<<"\n";
       }
       return;
  }
  return;
}
//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::PresFromRhoEg(Real rho, Real egas)
//  \brief Return gas pressure
Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {

  
  TempInvert(rho, LastTemp, egas, EOS::iE, EosData);
  LastTemp = EosData(EOS::iT);
 // if(rho<10.0){std::cout<<"Yes in presfromrhoeg  "<<rho<<"    "<<egas<<"   "<<LastTemp<<"\n";}
  return EosData(EOS::iP);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::EgasFromRhoP(Real rho, Real pres)
//  \brief Return internal energy density
Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {

  TempInvert(rho, LastTemp, pres, EOS::iP, EosData);

  LastTemp = EosData(EOS::iT);
 // if(rho<10.0){std::cout<<"Yes in egasfromrhoP  "<<rho<<"    "<<pres<<"   "<<LastTemp<<"\n";}
  return EosData(EOS::iE);
}

//----------------------------------------------------------------------------------------
//! \fn Real EquationOfState::AsqFromRhoP(Real rho, Real pres)
//  \brief Return adiabatic sound speed squared
Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {

  TempInvert(rho, LastTemp, pres, EOS::iP, EosData);
  
  LastTemp = EosData(EOS::iT);
 // if(rho<10.0){std::cout<<"Yes in asq  "<<rho<<"    "<<pres<<"   "<<LastTemp<<"\n";}
  return EosData(EOS::iAsq);
}

void EquationOfState::SevenFromRhoT(Real rho, Real T, AthenaArray<Real> &out) {
  QWData(rho, T, out);
}

Real EquationOfState::TFromRhoP(Real rho, Real pres) {

  TempInvert(rho, LastTemp, pres, EOS::iP, EosData);
  LastTemp = EosData(EOS::iT);
 // if(rho<10.0){std::cout<<"Yes in TfromrhoP  "<<rho<<"    "<<pres<<"   "<<LastTemp<<"\n";}
  return LastTemp;
}

Real EquationOfState::TFromRhoEgas(Real rho, Real egas) {

  TempInvert(rho, LastTemp, egas, EOS::iE, EosData);
  LastTemp = EosData(EOS::iT);
 // if(rho<10.0){std::cout<<"Yes in Tfromrhoegas  "<<rho<<"    "<<egas<<"   "<<LastTemp<<"\n";}
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
  T_floor = pin->GetOrAddReal("hydro", "T_floor", Tmin);
  T_ceil = pin->GetOrAddReal("hydro", "T_ceil", Tmax);
  LastTemp = pin->GetOrAddReal("hydro", "T0", std::sqrt(T_floor * T_ceil));
  nmax = pin->GetOrAddInteger("hydro", "nmax", 10000);
//  std::cout<<Ye<<" Ye\n";
  #ifdef MYDEBUG
    Real rho, temp,press;
    std::cout << "Input fluid parameters and retrieve EOS parameters." << '\n'
              << "Non-positive inputs will exit loop." << '\n';
    std::cout << "Input density (mass/volume): ";
    std::cin >> rho;
    std::cout << "Input temperature (K): ";
    std::cin >> temp;

    while (rho > 0 && std::isfinite(rho) && temp >0 && std::isfinite(temp)) {
      printf("d, t: %.16e, %.16e\n", rho, temp);
      QWData(rho, temp, EosData);
      Real p = EosData(EOS::iP);
      Real e = EosData(EOS::iE);
      Real a2 = EosData(EOS::iAsq);
      Real Te = TFromRhoEgas(rho, e);
      Real Tp = TFromRhoP(rho, p);
      printf("e(d, e)          , p(d, e)         , ASq(d, P)\n");
      printf("%.16e, %.16e, %.16e\n", e, p, a2);
      printf("T(d, e)/T-1, T(d, p)/T-1\n");
      printf("%.16e, %.16e\n", Te/temp - 1.0, Tp/temp - 1.0);
      printf("Root find %.16e \n",root_find(rho,p));
      std::cout << "Input density (mass/volume): ";
      std::cin >> rho;
      std::cout << "Input temperature (K): ";
      std::cin >> temp;
    }
    std::cout << std::endl;
  #endif
  return;
}
