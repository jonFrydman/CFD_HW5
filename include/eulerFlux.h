#ifndef EULERFLUX_H
#define EULERFLUX_H

#include "cellState.h"
#include "grid.h"
#include <vector>

vector<double> GenericFlux(grid &grd,  vector<cellState> &stencil);

vector<double> EastFlux(grid &grd, vector<cellState> &stencil);
vector<double> WestFlux(grid &grd, vector<cellState> &stencil);
vector<double> NorthFlux(grid &grd, vector<cellState> &stencil);
vector<double> SouthFlux(grid &grd, vector<cellState> &stencil);

#endif // EULERFLUX_H
