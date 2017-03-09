#ifndef NUMERICALFLUX_H_INCLUDED
#define NUMERICALFLUX_H_INCLUDED
// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.

#include "cellState.h"
#include "grid.h"
#include <vector>

using namespace std;
vector<double> GenericFlux(grid &grd, vector< vector<cellState> > &cellset, int i, int j, int delta_i, int delta_j);
vector<double> EastFlux(grid &grd, vector< vector<cellState> > &cellset, int i, int j);
vector<double> WestFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> NorthFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> SouthFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> JamesonViscocity(grid &grd, vector< vector<cellState> > &cellset,int i, int j);

//void IFT_IICD(int index, double dX, double dT, std::vector<double> & uT, std::vector<double> & uTplus);
//void IBT_2CD_I(int Xelements, std::vector<double>& CoefficientMinus, std::vector<double>& CoefficientPlus, std::vector<double>& Coefficient, std::vector<double>& Uplus, std::vector<double>& U);
// This is the end of the header guard
#endif // NUMERICALFLUX_H_INCLUDED
