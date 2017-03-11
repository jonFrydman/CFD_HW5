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

double nu_max(grid &grd, vector< vector<cellState> > &cellset,int i, int j, int delta_i, int delta_j);
vector<double> GenericJamesonViscocity(grid &grd, vector< vector<cellState> > &cellset,int i, int j, int delta_i, int delta_j);
vector<double> NorthJamVisc(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> SouthJamVisc(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> EastJamVisc(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> WestJamVisc(grid &grd, vector< vector<cellState> > &cellset,int i, int j);

vector<double> NorthFlux_AV(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> SouthFlux_AV(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> EastFlux_AV(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> WestFlux_AV(grid &grd, vector< vector<cellState> > &cellset,int i, int j);

vector<double> Residuals(grid &grd, vector< vector<cellState> > &cellset,int i, int j, double CFL);
vector<double> AlphaRK();
double Tau(grid &grd, vector< vector<cellState> > &cellset, int i, int j, double CFL);

vector< vector<cellState> > RK4(grid &grd, vector< vector<cellState> > &cellset, double CFL);
#endif // NUMERICALFLUX_H_INCLUDED
