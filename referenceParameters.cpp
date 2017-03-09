#include "referenceParameters.h"
//double rho_ref, P_ref, rhou_ref, rhov_ref, rhoE_ref, T_ref, H_ref, c_ref; //reference variables
//double gamma, cp, cv, R; //thermo constants

//Air Parameters
double R=8314./28.966; // air gas constant J/Kg-K
cp=1005.; // J/Kg-K
cv=cp-R;
gamma=cp/cv;
//reference variable (STATE 6)
rho_ref=1.293; // Kg/m^3
P_ref=101325;// Pa
ru_ref=sqrt(P_ref*rho_ref); // Kg/m^2-s volumetric-momentum
T_ref=P_ref/rho_ref; // K
rE_ref=P_ref/(gamma-1)+pow(ru_ref/rho_ref,2)/2; //internal energy KJ/m^3,
H_ref=(rE_ref+P_ref)/rho_ref;

c_ref=sqrt(gamma*P_ref/rho_ref);
