#include "cellState.h"
#include "grid.h"
#include <vector>
#include <iostream>
#include "eulerFlux.h"


vector<double> GenericFlux(grid &grd, vector< vector<cellState> > &cellset, int i, int j, int delta_i, int delta_j){
    double FSTAR1, FSTAR2, FSTAR3, FSTAR4;
    double GSTAR1, GSTAR2, GSTAR3, GSTAR4;

    FSTAR1 = 0.5*(cellset[i][j].F1() + cellset[i + delta_i][j+delta_j].F1());
	FSTAR2 = 0.5*(cellset[i][j].F2() + cellset[i + delta_i][j+delta_j].F2());
	FSTAR3 = 0.5*(cellset[i][j].F3() + cellset[i + delta_i][j+delta_j].F3());
	FSTAR4 = 0.5*(cellset[i][j].F4() + cellset[i + delta_i][j+delta_j].F4());
	GSTAR1 = 0.5*(cellset[i][j].G1() + cellset[i + delta_i][j+delta_j].G1());
	GSTAR2 = 0.5*(cellset[i][j].G2() + cellset[i + delta_i][j+delta_j].G2());
	GSTAR3 = 0.5*(cellset[i][j].G3() + cellset[i + delta_i][j+delta_j].G3());
	GSTAR4 = 0.5*(cellset[i][j].G4() + cellset[i + delta_i][j+delta_j].G4());

	std::vector<double> GENERIC_FLUX(4, 0.0);

	GENERIC_FLUX[0] = FSTAR1*grd.xWnorm[i+delta_i][j+delta_j] + GSTAR1*grd.yWnorm[i+delta_i][j+delta_j];
	GENERIC_FLUX[1] = FSTAR2*grd.xWnorm[i+delta_i][j+delta_j] + GSTAR2*grd.yWnorm[i+delta_i][j+delta_j];
	GENERIC_FLUX[2] = FSTAR3*grd.xWnorm[i+delta_i][j+delta_j] + GSTAR3*grd.yWnorm[i+delta_i][j+delta_j];
	GENERIC_FLUX[3] = FSTAR4*grd.xWnorm[i+delta_i][j+delta_j] + GSTAR4*grd.yWnorm[i+delta_i][j+delta_j];

	return GENERIC_FLUX;
}

vector<double> EastFlux(grid &grd, vector< vector<cellState> > &cellset, int i, int j){
    return GenericFlux(grd, cellset, i, j, 1, 0);
}
vector<double> WestFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j){
    return GenericFlux(grd, cellset, i, j, -1, 0);
}
vector<double> NorthFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j){
    return GenericFlux(grd, cellset, i, j, 0, 1);
}
vector<double> SouthFlux(grid &grd, vector< vector<cellState> > &cellset,int i, int j){
    return GenericFlux(grd, cellset, i, j, 0, -1);
}

