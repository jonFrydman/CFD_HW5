#ifndef NUMERICALFLUX_H_INCLUDED
#define NUMERICALFLUX_H_INCLUDED
// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.

#include "cellState.h"
#include "grid.h"
#include "eulerFlux.h"
#include <vector>

using namespace std;

double nu_max(vector<cellState> stencil);

vector<double> GenericJamesonViscocity(grid &grd, vector<cellState> &stencil);
vector<double> EastJamVisc(grid &grd,  vector<cellState> &stencil);
vector<double> WestJamVisc(grid &grd,  vector<cellState> &stencil);
vector<double> NorthJamVisc(grid &grd,  vector<cellState> &stencil);
vector<double> SouthJamVisc(grid &grd,  vector<cellState> &stencil);

vector<double> EastFlux_AV(grid &grd,  vector<cellState> &stencil);
vector<double> WestFlux_AV(grid &grd,  vector<cellState> &stencil);
vector<double> NorthFlux_AV(grid &grd,  vector<cellState> &stencil);
vector<double> SouthFlux_AV(grid &grd,  vector<cellState> &stencil);

vector<double> AirfoilFlux(grid & grd, vector< vector<cellState> > &cellset, int i);
vector<double> InletOutletFlux(grid &grd, vector< vector<cellState> > &cellset, int i);

#endif // NUMERICALFLUX_H_INCLUDED
