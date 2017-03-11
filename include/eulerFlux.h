#ifndef EULERFLUX_H
#define EULERFLUX_H

#include "cellState.h"
#include "grid.h"
#include <vector>

using namespace std;
vector<double> GenericFlux(grid &grd, vector< vector<cellState> > &cellset, int i, int j, int delta_i, int delta_j);
vector<double> EastFlux(grid &grd, vector< vector<cellState> > &cellset, int i, int j);
vector<double> WestFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> NorthFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j);
vector<double> SouthFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j);


#endif // EULERFLUX_H
