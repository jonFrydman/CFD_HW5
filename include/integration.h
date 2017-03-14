#ifndef INTEGRATION_H_INCLUDED
#define INTEGRATION_H_INCLUDED

#include<vector>
#include "grid.h"
#include "cellState.h"
#include "numericalFlux.h"
#include "eulerFlux.h"

using namespace std;

vector<double> Residuals(grid &grd,  vector< vector<cellState> > &cellset, int i, int j);
double Tau(grid &grd, cellState cell, double CFL);
vector<double> AlphaRK();
vector< vector<cellState> > RK4(grid &grd, vector< vector<cellState> > &cellset, double CFL);

vector<cellState> stencilEW(grid & grd, vector< vector<cellState> > & cellset, int i, int j);
vector<cellState> stencilNS(grid & grd, vector< vector<cellState> > & cellset, int i, int j);

#endif // INTEGRATION_H_INCLUDED
