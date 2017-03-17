#ifndef REFERENCEPARAMETERS_H_INCLUDED
#define REFERENCEPARAMETERS_H_INCLUDED

#include<math.h>
//Problem Parameters
double M_ref=0.85;
double CFL=2.4; // p578  or hrich v1 says less than CFL=2.8

//Air Properties at STP
double R=8314./28.966; // air gas constant J/Kg-K
double cp=1005.; // J/Kg-K
double cv=cp-R;
double gamma=cp/cv;

//Reference pressure and density at STP (NB: stagnation pressures at different values);
double rho_ref=1.293; // Kg/m^3
double P_ref=101325;// Pa
double T_ref=P_ref/rho_ref; // K
double c_ref=sqrt(gamma*P_ref/rho_ref);
double s_ref=cv*log(T_ref)-R*log(rho_ref);

double rhoE_ref=P_ref/(gamma-1)+0.5*rho_ref*pow(M_ref*c_ref,2); //internal energy KJ/m^3,
double H_ref=(rhoE_ref+P_ref)/rho_ref;
double speed_ref=M_ref*c_ref;

double alpha2=1.0/4;
double alpha4=1.0/256;

#endif // REFERENCEPARAMETERS_H_INCLUDED
